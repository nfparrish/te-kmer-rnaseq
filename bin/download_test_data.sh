#!/usr/bin/env bash
set -euo pipefail
#
# Download a subset of E-MTAB-3788 (LINE-1 knockdown RNA-seq) from ENA.
# ENA project: ERP011233
#
# We download MCF-7 samples only (4 conditions × ~3 replicates) to keep
# the validation manageable.  Full sample list is in the SDRF:
# https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-3788
#
# Usage:  bash bin/download_test_data.sh [output_dir]

OUTDIR="${1:-test_data}"
mkdir -p "$OUTDIR"

# Representative MCF-7 run accessions from ERP011233:
#   untreated:   ERR1018546, ERR1018547
#   sh960 (ctrl): ERR1018548, ERR1018549, ERR1018550
#   sh1083 (KD):  ERR1018551, ERR1018552, ERR1018553
#   sh1085 (KD):  ERR1018554, ERR1018555, ERR1018556
RUNS=(
    ERR1018546 ERR1018547
    ERR1018548 ERR1018549 ERR1018550
    ERR1018551 ERR1018552 ERR1018553
    ERR1018554 ERR1018555 ERR1018556
)

echo "Downloading ${#RUNS[@]} paired-end run(s) into ${OUTDIR}/ ..."

for RUN in "${RUNS[@]}"; do
    PREFIX="${RUN:0:6}"
    SUBDIR="${RUN}"
    URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${PREFIX}/${SUBDIR}"

    echo "  → ${RUN}"
    wget -q -P "$OUTDIR" "${URL}/${RUN}_1.fastq.gz" || \
        echo "    WARN: could not download ${RUN}_1.fastq.gz – check ENA"
    wget -q -P "$OUTDIR" "${URL}/${RUN}_2.fastq.gz" || \
        echo "    WARN: could not download ${RUN}_2.fastq.gz – check ENA"
done

echo "Done.  Point --reads '${OUTDIR}/*_{1,2}.fastq.gz' at the pipeline."
