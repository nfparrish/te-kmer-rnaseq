#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * te-kmer-rnaseq: Quantify transposable element expression via k-mer
 * decomposition of reads mapped to Dfam TE family consensus sequences.
 *
 * Conceptually inspired by Kojima et al. 2021 hervk_kmers pipeline
 * (WGS k-mer counting for HERV-K), adapted here for RNA-seq and
 * generalised to all TE families with k-means clustering.
 */

// ── Parameter defaults (override in nextflow.config or CLI) ────────
params.reads        = null          // glob to paired-end FASTQ, e.g. 'data/*_{1,2}.fastq.gz'
params.species      = 'human'       // 'human' or 'mouse'
params.k            = 50            // k-mer length
params.max_mismatch = 10            // max mismatches tolerated during mapping
params.max_clip     = 10            // max soft-clip bases tolerated
params.n_clusters   = 50            // k-means clusters
params.te_consensus = null          // optional pre-built TE consensus FASTA
params.outdir       = 'results'
params.dfam_url     = 'https://www.dfam.org/releases/current/families/Dfam_curatedonly.h5.gz'

// ── Validate ───────────────────────────────────────────────────────
if (!params.reads) { error "Please supply --reads 'path/to/*_{1,2}.fastq.gz'" }

// ── Channel: paired-end reads ──────────────────────────────────────
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .set { ch_reads }

// ── Process: Obtain TE consensus FASTA from Dfam ───────────────────
process FETCH_DFAM_CONSENSUS {
    tag "${params.species}"
    publishDir "${params.outdir}/reference", mode: 'copy'

    output:
    path "te_consensus.fa", emit: fasta

    script:
    if (params.te_consensus)
        """
        cp ${params.te_consensus} te_consensus.fa
        """
    else
        // Use famdb.py from Dfam Tools to export consensus for the target species.
        // Alternatively, download the curated FASTA release directly.
        """
        #!/usr/bin/env bash
        set -euo pipefail

        SPECIES="${params.species == 'human' ? 'Homo sapiens' : 'Mus musculus'}"

        # Try Dfam curated FASTA (lighter than HDF5).
        # Fall back to a direct download of curated consensi.
        if command -v famdb.py &>/dev/null; then
            famdb.py -i Dfam.h5 families \
                --format fasta_name \
                --include-class-in-name \
                --ancestors --descendants "\$SPECIES" \
                > te_consensus.fa
        else
            # Download curated consensi FASTA from Dfam FTP
            curl -fSL "https://www.dfam.org/releases/current/families/Dfam_curatedonly.fasta.gz" \
                -o dfam.fa.gz
            gunzip dfam.fa.gz

            # Filter to target species (Dfam FASTA headers contain OS: tag)
            python3 - <<'PYEOF'
import re, sys
keep = False
target = "${params.species == 'human' ? 'Homo sapiens' : 'Mus musculus'}"
with open("dfam.fa") as fin, open("te_consensus.fa", "w") as fout:
    for line in fin:
        if line.startswith(">"):
            keep = target.lower() in line.lower()
        if keep:
            fout.write(line)
PYEOF
            # Ensure non-empty
            if [ ! -s te_consensus.fa ]; then
                echo "WARN: species filter yielded empty file; using full Dfam set" >&2
                mv dfam.fa te_consensus.fa 2>/dev/null || cp /dev/null te_consensus.fa
            fi
        fi
        """
}

// ── Process: Build bowtie2 index ───────────────────────────────────
process BUILD_BOWTIE2_INDEX {
    tag "bt2_index"

    input:
    path fasta

    output:
    path "te_idx*", emit: index

    script:
    """
    bowtie2-build --threads ${task.cpus} ${fasta} te_idx
    """
}

// ── Process: Read QC & trimming ────────────────────────────────────
process FASTP {
    tag "${sample_id}"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: '*.json'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_{1,2}.fastq.gz"), emit: trimmed
    tuple val(sample_id), path("${sample_id}_fastp.json"),              emit: json

    script:
    """
    fastp \
        -i ${reads[0]} -I ${reads[1]} \
        -o ${sample_id}_trimmed_1.fastq.gz \
        -O ${sample_id}_trimmed_2.fastq.gz \
        -j ${sample_id}_fastp.json \
        --thread ${task.cpus} \
        --qualified_quality_phred 20 \
        --length_required 50
    """
}

// ── Process: Count total reads for depth normalisation ─────────────
process COUNT_TOTAL_READS {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(json)

    output:
    tuple val(sample_id), env(TOTAL_READS), emit: depth

    script:
    """
    TOTAL_READS=\$(python3 -c "
import json, sys
with open('${json}') as f:
    d = json.load(f)
print(d['summary']['after_filtering']['total_reads'])
")
    """
}

// ── Process: Map reads to TE consensus (error-tolerant) ────────────
process MAP_TO_TE_CONSENSUS {
    tag "${sample_id}"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}.te.bam"), path("${sample_id}.te.bam.bai"), emit: bam

    script:
    // --very-sensitive-local: high sensitivity, local alignment handles partial TE matches
    // --no-mixed: both mates must align for concordant
    // -N 1 --np 0: tolerate seed mismatches, no penalty for Ns
    // --score-min: relax threshold proportional to read length
    """
    bowtie2 \
        --very-sensitive-local \
        -N 1 --np 0 \
        --score-min L,0,0.5 \
        --no-unal \
        -p ${task.cpus} \
        -x te_idx \
        -1 ${reads[0]} -2 ${reads[1]} \
    | samtools view -@ ${task.cpus} -bS -F 4 - \
    | samtools sort -@ ${task.cpus} -o ${sample_id}.te.bam

    samtools index ${sample_id}.te.bam
    """
}

// ── Process: Extract k-mers from mapped reads ──────────────────────
process EXTRACT_KMERS {
    tag "${sample_id}"
    publishDir "${params.outdir}/kmers", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_kmer_counts.tsv.gz"), emit: counts

    script:
    """
    extract_kmers.py \
        --bam ${bam} \
        --k ${params.k} \
        --max-mismatch ${params.max_mismatch} \
        --max-clip ${params.max_clip} \
        --output ${sample_id}_kmer_counts.tsv.gz
    """
}

// ── Process: Merge k-mer counts across all samples ─────────────────
process MERGE_KMER_COUNTS {
    tag "merge"
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    path kmer_files      // collected list
    path depth_file      // TSV of sample_id\ttotal_reads

    output:
    path "kmer_matrix_raw.tsv.gz",  emit: raw_matrix
    path "kmer_matrix_cpm.tsv.gz",  emit: cpm_matrix

    script:
    """
    normalize_and_cluster.py merge \
        --kmer-files ${kmer_files} \
        --depth-file ${depth_file} \
        --out-raw kmer_matrix_raw.tsv.gz \
        --out-cpm kmer_matrix_cpm.tsv.gz
    """
}

// ── Process: K-means clustering of k-mers ──────────────────────────
process CLUSTER_KMERS {
    tag "cluster"
    publishDir "${params.outdir}/clusters", mode: 'copy'

    input:
    path cpm_matrix

    output:
    path "kmer_clusters.tsv",              emit: clusters
    path "cluster_expression.tsv",         emit: expression
    path "cluster_elbow.png",              emit: elbow_plot
    path "cluster_pca.png",                emit: pca_plot

    script:
    """
    normalize_and_cluster.py cluster \
        --cpm-matrix ${cpm_matrix} \
        --n-clusters ${params.n_clusters} \
        --out-clusters kmer_clusters.tsv \
        --out-expression cluster_expression.tsv \
        --out-elbow cluster_elbow.png \
        --out-pca cluster_pca.png
    """
}

// ── Process: Collect depth info into a single file ─────────────────
process COLLECT_DEPTHS {
    tag "depths"

    input:
    val depths   // list of [sample_id, total_reads]

    output:
    path "sample_depths.tsv", emit: depth_file

    script:
    def lines = depths.collect { "${it[0]}\t${it[1]}" }.join("\\n")
    """
    printf "sample_id\\ttotal_reads\\n${lines}\\n" > sample_depths.tsv
    """
}

// ── Workflow ────────────────────────────────────────────────────────
workflow {

    // 1. Obtain / prepare TE consensus FASTA
    FETCH_DFAM_CONSENSUS()

    // 2. Build bowtie2 index
    BUILD_BOWTIE2_INDEX(FETCH_DFAM_CONSENSUS.out.fasta)

    // 3. Trim reads
    FASTP(ch_reads)

    // 4. Count total reads for normalisation
    COUNT_TOTAL_READS(FASTP.out.json)

    // 5. Map trimmed reads to TE consensus sequences
    MAP_TO_TE_CONSENSUS(FASTP.out.trimmed, BUILD_BOWTIE2_INDEX.out.index.collect())

    // 6. Extract k-mers from mapped reads
    EXTRACT_KMERS(MAP_TO_TE_CONSENSUS.out.bam)

    // 7. Collect depths
    COLLECT_DEPTHS(COUNT_TOTAL_READS.out.depth.collect())

    // 8. Merge k-mer counts across samples
    MERGE_KMER_COUNTS(
        EXTRACT_KMERS.out.counts.map { it[1] }.collect(),
        COLLECT_DEPTHS.out.depth_file
    )

    // 9. Cluster k-mers
    CLUSTER_KMERS(MERGE_KMER_COUNTS.out.cpm_matrix)
}
