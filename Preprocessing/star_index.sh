#!/bin/sh -x

# Set name of the job
#SBATCH --job-name=Index_Customref    # Job name

# ── Paths ─────────────────────────────────────────────────────────────────────
directory="PATH/TO/genome_reference_directory"

genome_fasta="${directory}/genome.primary_assembly.fa"
annotation_gtf="${directory}/annotation.primary_assembly.gtf"
# ──────────────────────────────────────────────────────────────────────────────

# ── Build STAR genome index ───────────────────────────────────────────────────
mkdir -p "${directory}/reference"

STAR \
    --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir "${directory}/reference" \
    --genomeFastaFiles "${genome_fasta}" \
    --sjdbGTFfile "${annotation_gtf}"
