#!/bin/bash	
# Expansion
find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "java -jar -Xmx8g /usr/local/genome/picard_latest/picard.jar FastqToSam  FASTQ={}.fq   O=/40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}_unmapped.bam SM={}" 
find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "/usr/local/genome/bin/STAR-2.7.1a --runMode alignReads --genomeDir /40/Cruchaga_Data/bulkRNASeq/References/GRCh38/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v33_oh100 --twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 --readFilesIn /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}_unmapped.bam --readFilesType SAM SE --readFilesCommand samtools view -h --outFileNamePrefix /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}. --outSAMstrandField intronMotif  --outFilterIntronMotifs None --alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted --outSAMunmapped Within --genomeLoad NoSharedMemory --limitGenomeGenerateRAM 31000000000 --genomeChrBinNbits 14 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimOutType WithinBAM SoftClip --chimMainSegmentMultNmax 1 --outSAMattributes NH HI AS nM NM ch"
find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "samtools sort --threads 2 -o /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.Aligned.sortedByCoord.out.bam /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.Aligned.out.bam"
find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "samtools index /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.Aligned.sortedByCoord.out.bam"
find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "java -jar -Xmx8g /usr/local/genome/picard-2.26.0/picard.jar CollectRnaSeqMetrics I= /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.Aligned.sortedByCoord.out.bam O= /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.RNA_Metrics.txt REF_FLAT= /40/Cruchaga_Data/bulkRNASeq/References/GRCh38/GRCh38.refFlat.txt RIBOSOMAL_INTERVALS= /40/Cruchaga_Data/bulkRNASeq/References/GRCh38/gencode.v33.rRNA.interval_list.txt STRAND_SPECIFICITY=NONE TMP_DIR= /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/tmp MAX_RECORDS_IN_RAM=150000"
find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "java -jar -Xmx8g /usr/local/genome/picard-2.26.0/picard.jar CollectAlignmentSummaryMetrics I=/40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.Aligned.sortedByCoord.out.bam ADAPTER_SEQUENCE=null O=/40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}_Summary_metrics.txt TMP_DIR= /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/tmp"
find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "/usr/local/genome/salmon-0.11.3/bin/salmon quant --threads 2 -t /40/Cruchaga_Data/bulkRNASeq/References/GRCh38/gencode.v33.transcripts.ERCC_4_salmon.fa -a /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.Aligned.toTranscriptome.out.bam --libType A --seqBias --gcBias --geneMap /40/Cruchaga_Data/bulkRNASeq/References/GRCh38/gencode.v33.primary_assembly.annotation.gtf -o /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.salmonp"
done


find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "java -jar -Xmx8g /usr/local/genome/picard-2.26.0/picard.jar CollectRnaSeqMetrics I= /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.Aligned.sortedByCoord.out.bam O= /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.RNA_Metrics.txt REF_FLAT= /40/Cruchaga_Data/bulkRNASeq/References/GRCh38/GRCh38.refFlat.txt RIBOSOMAL_INTERVALS= /40/Cruchaga_Data/bulkRNASeq/References/GRCh38/gencode.v33.rRNA.interval_list.txt STRAND_SPECIFICITY=NONE TMP_DIR= /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/tmp MAX_RECORDS_IN_RAM=150000"

find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 "java -jar -Xmx8g /usr/local/genome/picard-2.26.0/picard.jar CollectAlignmentSummaryMetrics I=/40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}.Aligned.sortedByCoord.out.bam ADAPTER_SEQUENCE=null O=/40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/{}_Summary_metrics.txt TMP_DIR= /40/AD/CellFree_Seq_Data/01-cfRNA/Expansion_hg38/tmp"

find *.fq |  rev | cut -c 4- | rev | uniq | parallel -j 6 

# HH K99 DoD

mkdir -p $out_dir/$sample
#sample_full_name="${1}.star.bam"


java -jar -Xmx40g /opt/picard-tools/picard.jar FastqToSam FASTQ=$raw_data/$sample.fastq.gz O=$out_dir/$sample/$sample_unmapped.bam SM=$sample


/opt/STAR-2.7.1a/bin/Linux_x86_64/STAR --runMode alignReads --genomeDir $path2Ref_HHC/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v33_oh100 \
--twopassMode Basic --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType \
BySJout --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --limitSjdbInsertNsj 1200000 \
--readFilesIn $out_dir/$sample/$sample_unmapped.bam --readFilesType SAM SE --readFilesCommand samtools view -h \
--outFileNamePrefix $out_dir/$sample/$sample. --outSAMstrandField intronMotif --outFilterIntronMotifs None \
--alignSoftClipAtReferenceEnds Yes --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted \
--outSAMunmapped Within --genomeLoad NoSharedMemory --limitGenomeGenerateRAM 31000000000 --genomeChrBinNbits 14 \
--chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimOutType WithinBAM SoftClip --chimMainSegmentMultNmax 1 \
--outSAMattributes NH HI AS nM NM ch

samtools sort --threads 12 -o $out_dir/$sample/$sample.Aligned.sortedByCoord.out.bam $out_dir/$sample/$sample.Aligned.out.bam

samtools index $out_dir/$sample/$sample.Aligned.sortedByCoord.out.bam

java -jar -Xmx40g /opt/picard-tools/picard.jar CollectRnaSeqMetrics I=$out_dir/$sample/$sample.Aligned.sortedByCoord.out.bam \
O=$out_dir/$sample/$sample.RNA_Metrics.txt REF_FLAT=$path2Ref/GRCh38.refFlat.txt \
RIBOSOMAL_INTERVALS=$path2Ref/gencode.v33.rRNA.interval_list.txt \
STRAND_SPECIFICITY=NONE TMP_DIR=./tmp MAX_RECORDS_IN_RAM=150000

java -jar -Xmx40g /opt/picard-tools/picard.jar CollectAlignmentSummaryMetrics I=$out_dir/$sample/$sample.Aligned.sortedByCoord.out.bam \
ADAPTER_SEQUENCE=null O=$out_dir/$sample/$sample_Summary_metrics.txt TMP_DIR=./tmp

java -jar -Xmx40g /opt/picard-tools/picard.jar MarkDuplicates I=$out_dir/$sample/$sample.Aligned.sortedByCoord.out.bam \
O=$out_dir/$sample/$sample.Aligned.sortedByCoord.out.md.bam PROGRAM_RECORD_ID=null M=$out_dir/$sample/$sample.marked_dp_metrics.txt \
TMP_DIR=./tmp ASSUME_SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=150000 OPTICAL_DUPLICATE_PIXEL_DISTANCE=100

/opt/salmon-latest_linux_x86_64/bin/salmon quant --threads 12 -t $path2Ref/gencode.v33.transcripts.ERCC_4_salmon.fa \
-a $out_dir/$sample/$sample.Aligned.toTranscriptome.out.bam --libType A --seqBias --gcBias --geneMap $path2Ref/gencode.v33.GRCh38.annotation.ERCC92.gtf \
-o $out_dir/$sample/$sample.salmon
