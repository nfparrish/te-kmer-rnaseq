# te-kmer-rnaseq

Nextflow pipeline that quantifies **transposable element (TE) expression heterogeneity** from RNA-seq by decomposing TE-mapped reads into k-mers and clustering them.

## Approach

1. **Map** paired-end RNA-seq reads to Dfam TE family consensus sequences (human/mouse) with bowtie2 in error-tolerant local-alignment mode.
2. **Decompose** mapped reads (and their mates) into sliding-window k-mers (default k = 50), filtering by mismatch/clip thresholds (inspired by [Kojima et al. 2021](https://github.com/GenomeImmunobiology/Kojima_et_al_2021_hervk_kmers)).
3. **Normalise** k-mer counts to CPM using total sequencing depth.
4. **Cluster** k-mers (k-means) by normalised expression profile — k-mers in the same cluster likely reflect the same TE transcript sequence (defined by sequence, not locus).
5. **Output** cluster-level expression for downstream association with genetic variation.

Because TEs are repetitive and polymorphic across individuals, this reference-free k-mer strategy captures TE sequence-specific expression that conventional genome-mapping pipelines miss.

## Quick start

```bash
# Install: nextflow, bowtie2, samtools, fastp, python3 (pysam, pandas, scikit-learn, matplotlib)

nextflow run main.nf \
    --reads 'data/*_{1,2}.fastq.gz' \
    --species human \
    --k 50 \
    --n_clusters 50 \
    --outdir results
```

## Validation with E-MTAB-3788

The [E-MTAB-3788](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-3788) dataset contains 2×150 bp paired-end RNA-seq from MCF-7 and 2102Ep cells with LINE-1 (L1 ORF1) shRNA knockdown. LINE-1 k-mer clusters should show reduced expression in knockdown vs control.

```bash
bash bin/download_test_data.sh   # downloads MCF-7 runs from ENA
nextflow run main.nf -profile test
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | — | Glob to paired-end FASTQ files |
| `--species` | `human` | `human` or `mouse` |
| `--k` | `50` | K-mer length |
| `--max_mismatch` | `10` | Max mismatches per read |
| `--max_clip` | `10` | Max soft-clipped bases |
| `--n_clusters` | `50` | Number of k-means clusters |
| `--te_consensus` | — | Pre-built TE consensus FASTA (skips Dfam download) |

## Outputs

- `results/kmers/` — per-sample k-mer count tables
- `results/merged/` — raw and CPM-normalised k-mer × sample matrices
- `results/clusters/` — k-mer cluster assignments, per-cluster expression, elbow plot, PCA plot

## License

MIT
