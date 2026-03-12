#!/bin/sh -x

# Set name of the job
#SBATCH --job-name=Map_Mm_CDKL5    # Job name

# Activate conda environment
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate YOUR_CONDA_ENV   # Your conda environment name

# ── Paths ─────────────────────────────────────────────────────────────────────
directory="PATH/TO/project_directory"
data_directory="PATH/TO/project_directory/fastq"

bbduk="PATH/TO/bbmap/bbduk.sh"
bbmap_directory="PATH/TO/bbmap/resources"

STAR_indices="PATH/TO/STAR/genome_index"

annotation_gtf="PATH/TO/annotation.gtf"
# ──────────────────────────────────────────────────────────────────────────────

# ── Trimming ──────────────────────────────────────────────────────────────────
mkdir -p "${directory}/Trimmed"

for sample in "${data_directory}"/*.fastq.gz; do
    last_part=$(basename "$sample")
    $bbduk \
        in="$sample" \
        out="${directory}/Trimmed/${last_part}" \
        ref="${bbmap_directory}/polyA.fa.gz,${bbmap_directory}/truseq.fa.gz" \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20
done

# ── STAR Alignment ────────────────────────────────────────────────────────────
mkdir -p "${directory}/STARoutput"

for clean_sample in "${directory}/Trimmed"/*.fastq.gz; do
    last_part=$(basename "$clean_sample")
    STAR --runThreadN 8 \
         --genomeDir "${STAR_indices}" \
         --readFilesIn "${clean_sample}" \
         --outFilterType BySJout \
         --outFilterMultimapNmax 20 \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverLmax 0.1 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --outSAMattributes NH HI NM MD \
         --readFilesCommand "gunzip -c" \
         --outReadsUnmapped Fastx \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix "${directory}/STARoutput/${last_part%%_clean*}_"
done

# ── SAMtools Indexing ─────────────────────────────────────────────────────────
for bamfile in "${directory}/STARoutput"/*.out.bam; do
    samtools index "${bamfile}"
done

# ── Read Counting (HTSeq) ─────────────────────────────────────────────────────
# Note: stranded=yes for DGE reads
mkdir -p "${directory}/ReadCounts"

for bamfile in "${directory}/STARoutput"/*.out.bam; do
    sample_name=$(basename "$bamfile")
    htseq-count \
        -m intersection-nonempty \
        -s yes \
        -f bam \
        -r pos \
        "$bamfile" \
        "$annotation_gtf" \
        > "${directory}/ReadCounts/${sample_name}_Read_Counts.txt"
done
