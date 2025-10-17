THREADS=8
GENOME_DIR=~/task_1/star_index/GRCh38.99.chr21
FA=~/task_1/ASimulatoR/inst/extdata/21.fa
GTF=~/task_1/ASimulatoR/inst/extdata/Homo_sapiens.GRCh38.99.21.gtf
FASTQ_DIR=~/task_1/ASimulatoR/myout/2025-10-12-13:49:36.810135_maxGenes0_SeqDepth2e+06_errRate0.001_readlen76_multiEventsPerExonTRUE_probsAsFreqFALSE
R1=$FASTQ_DIR/sample_01_1.fastq
R2=$FASTQ_DIR/sample_01_2.fastq
OUTDIR=~/task_1/star_align/sample_01
mkdir -p "$GENOME_DIR" "$OUTDIR"

STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$FA" --sjdbGTFfile "$GTF" --sjdbOverhang 75

STAR --runThreadN $THREADS --genomeDir "$GENOME_DIR" \
     --readFilesIn "$R1" "$R2" \
     --outFileNamePrefix "$OUTDIR/" \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts
