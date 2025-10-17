THREADS=16

GENOME_DIR=~/task_1/star_index/GRCh38.115.CHR
FA_DIR=~/task_1/data/ensembl_release_115
GTF=~/task_1/data/ensembl_release_115/Homo_sapiens.GRCh38.115.chr.gtf
FASTQ_DIR=~/task_1/ASimulatoR/big_out11/2025-10-16-17:06:34.316872_maxGenes0_SeqDepth2e+06_errRate0.001_readlen76_multiEventsPerExonTRUE_probsAsFreqFALSE
R1=$FASTQ_DIR/sample_01_1.fastq
R2=$FASTQ_DIR/sample_01_2.fastq
OUTDIR=~/task_1/star_align_big/sample_01
mkdir -p "$GENOME_DIR" "$OUTDIR"

STAR --runThreadN $THREADS \
    --runMode genomeGenerate \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles $FA_DIR/{1..22}.fa $FA_DIR/X.fa $FA_DIR/Y.fa \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang 75

STAR --runThreadN $THREADS \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$R1" "$R2" \
    --outFileNamePrefix "$OUTDIR/" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts
