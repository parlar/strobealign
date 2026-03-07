#!/bin/bash
set -euo pipefail

REF=/home/parlar_ai/dev/strobealign/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
R1=/home/parlar_ai/dev/sv_tools/data/giab_hg38/fastq/HG002.novaseq.pcr-free.40x.R1.fastq.gz
R2=/home/parlar_ai/dev/sv_tools/data/giab_hg38/fastq/HG002.novaseq.pcr-free.40x.R2.fastq.gz
OUTDIR=/home/parlar_ai/dev/sv_tools/data/giab_hg38/HG002
OUT=${OUTDIR}/HG002.GRCh38.novaseq-40x.strobealign-sv.two-pass.bam
LOG=${OUTDIR}/HG002.GRCh38.novaseq-40x.strobealign-sv.two-pass.log

/home/parlar_ai/dev/strobealign/build/strobealign \
  -t 20 \
  --two-pass \
  --eqx \
  --supp=3 \
  --sv-tags \
  --use-index \
  "$REF" \
  "$R1" \
  "$R2" \
  2>"$LOG" \
  | samtools sort -@4 -o "$OUT"

samtools index "$OUT"

echo "Done. BAM: $OUT  Log: $LOG"
