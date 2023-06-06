#!/bin/bash


############################################################################################################################################################################################
# HI THERE!
#
# Analyze bulk ATACseq data: you do NOT need to change any of this script, just define args
#
# To run this pipe with your own samples, please just provide the corresponding 2 arguments: ("$1") path to read1 and ("$2") path to output folder #######?¿??
# They have tag at the end: 	# Specificed by USER
#
# Example:  Run and send to background 
	# ./1.analyzeRLs_AllReads.sh /dd_rundata/novaseq/Runs/190608_A00690_H7V22DRXX_004/Unaligned/snATAC/RL1753_2019_06_08_bulkATAC_20190527_C4593_50000nuclei_37deg_5ulTn5_S1_R1_001.fastq.gz /home/aalvarez/Bulk_ATAC_QCmetrics_tool/output > RL1760.log  2>&1 &  
#
# Feed this script with all your RLs at once (loop) or one by one if you wanna paralelize the analysis
#
# You can set the same out dir for all your RLs, results are not going to overlap, the script will use the same folder tree
#
# NOTE: PARAMs are set for PairEnd data!!!
#
# very very important questions: aalvarezf@cnic.es
############################################################################################################################################################################################

# If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

###### Inputs: 

R1="$1" # Specificed by USER: arg.1
R2=${R1/_R1./_R2.}
RL=$( basename $R1 | sed 's/_.*//g')

echo -e "Working with RL: $RL\nR1: $R1\nR2: $R2"

###### Outputs: 

outD="$2" # Specificed by USER: arg.2

dir0=$outD"/0.LOGs/"
dir1=$outD"/1.QCs/"
dir2=$outD"/2.ProcessedReads"
dir3=$outD"/3.QCsRmAdaptors"
dir4=$outD"/4.Alignment"

#tool="$3" # Specificed by USER: arg.3
tool="/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool" # Specificed by USER: arg.3

mkdir -p $dir0 $dir1 $dir2 $dir3 $dir4 

###### Utilis: 
dataD=$tool"/data"  
FASTQC="/home/aalvarez/Tools/FastQC/fastqc"
#BBDUK="/home/sbuckberry/working_data_01/bin/bbmap/bbduk2.sh"
CUTADAPT="/home/dmakosa/working_data_01/apps/miniconda3/envs/10xmethylomes/bin/cutadapt"
sambamba="/home/dmakosa/working_data_01/apps/miniconda3/envs/10xmethylomes/bin/sambamba"
parallel="/home/sbuckberry/working_data_01/bin/parallel"
blacklist=$dataD"/ENCFF000KJP.bed"
chromSizes=$dataD"/hg19.chromSizes"
promoters=$dataD"/promoterhg19_2kb.bed" # adapted, 2KB long
chrBed=$dataD"/chromS.bed"
index="/home/sbuckberry/working_data_01/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
black_M=$dataD"/M_size.bed"


###### Aanalyze:

# 1.QCs raw reads
#echo -e "\t1. QCs"
#$FASTQC -o $dir1 $R1 $R2 > ${dir0}/logPreQCs_${RL}

# 2.Clean Adaptors	
#echo -e "\t2. Clean Adaptors"
#$BBDUK in=$R1 in2=$R2 out=$dir2"/"$RL"_R1_clean.fq.gz" out2=$dir2"/"$RL"_R2_clean.fq.gz" rliteral=GCGATCGAGGACGGCAGATGTGTATAAGAGACAG,CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG ktrim=r mink=3 threads=3 overwrite=true > ${dir0}"/logAdap"$RL 
#$CUTADAPT -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o $dir2"/"$RL"_R1_clean.fq.gz" -p $dir2"/"$RL"_R2_clean.fq.gz" $R1 $R2 > ${dir0}"/logAdap"$RL 
	
# 3. QCs after clean reads
#echo -e "\t3. QCs after adaptors"
#$FASTQC -o $dir3  $dir2"/"$RL"_R1_clean.fq.gz"  $dir2"/"$RL"_R2_clean.fq.gz" > ${dir0}/logQCs_${RL}

#4. Bowtie2 alignment::care or -T option, do not overlap temp files
#echo -e "\t4. Bowtie and samtools sort" 
#(bowtie2 -q --threads 4 -X2000 -x $index -1 $dir2"/"$RL"_R1_clean.fq.gz"  -2 $dir2"/"$RL"_R2_clean.fq.gz" | samtools view -bSu - | samtools sort -T ${RL}sorted - > $dir4"/"${RL}.bam) >  ${dir0}/logBowtie${RL}

# 5. STATS
#echo -e "\t5. Flag Stats"
#samtools flagstat $dir4"/"${RL}.bam  > $dir4"/"${RL}_flagStats

# 6. Count MT reads
#echo -e "\t6. Count MT reads"
#samtools index  $dir4"/"${RL}.bam # needed for -c
#samtools view -c $dir4"/"${RL}.bam chrM > $dir4"/"${RL}_mitoReads


# 7. Proper Pairs
#echo -e "\t7. Data is PE ---> $pairs\n\tFilter alignment, proper Pairs and stats" 
$sambamba view -t 4 -f bam -F "proper_pair" $dir4"/"${RL}.bam -o $dir4"/"${RL}_pairs.bam 
samtools flagstat $dir4"/"${RL}_pairs.bam > $dir4"/"${RL}_pairs_flagStats


# 8. PCR dups
echo -e "\t8. Remove PCR dups and stats"
$sambamba markdup -p -t 4 -r --hash-table-size=1000000  $dir4"/"${RL}_pairs.bam $dir4"/"${RL}_pairs_dedup.bam
samtools flagstat $dir4"/"${RL}_pairs_dedup.bam > $dir4"/"${RL}_pairs_dedup_flagStats


# 9. Remove blacklist 
echo -e "\t9. Intercept BlackListed regions"
samtools view -q 30 -L $blacklist -U $dir4"/"${RL}_pairs_dedup_filt.bam -b $dir4"/"${RL}_pairs_dedup.bam > $dir4"/"${RL}_pairs_dedup_crap
samtools flagstat $dir4"/"${RL}_pairs_dedup_filt.bam > $dir4"/"${RL}_pairs_dedup_filt_flagStats
	

# 10. Remove MT reads
samtools view -L $black_M -U $dir4"/"${RL}_pairs_dedup_filt_noMT.bam -b $dir4"/"${RL}_pairs_dedup_filt.bam > $dir4"/"${RL}_pairs_dedup_2crap
samtools flagstat $dir4"/"${RL}_pairs_dedup_filt_noMT.bam > $dir4"/"${RL}_pairs_dedup_filt_noMT_flagStats
	
# Coverage --> may break the pipeline!
samtools index $dir4"/"${RL}_pairs_dedup_filt_noMT.bam 
/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/bamCoverage --binSize 50 -p 5 -b $dir4"/"${RL}_pairs_dedup_filt_noMT.bam -o $dir4"/"${RL}_pairs_dedup_filt_noMT_coverage.bigwig # maybe problems here, depending on size

	# 11. RMs		
echo -e "\t11. Remove Intermediate Files, some bam also"
rm $dir4"/"${RL}_*crap	
rm $dir4"/"${RL}_pairs.bam*
rm $dir4"/"${RL}_pairs_dedup.bam*
rm $dir4"/"${RL}_pairs_dedup_filt.bam

# Just Keep original alignment + final filtered bam + index final + flagStats all + coverage.bw final
fBam=$dir4"/"${RL}_pairs_dedup_filt_noMT.bam

# 12. Reads in Promotors
echo -e "\t12. Counting reads in promoters"
samtools view -c -L $promoters -b $fBam > $dir4/$RL"_ReadsInProm"
	
# 13. MACS2 NFR peaks
echo -e "\t13. MACS2:::params for Tn5 insert sites, params ENCODE pipe for NFRs\n\tCalcule Reads in NFR peaks\n\tCalcule Reads in NFR peaks overlapping promotors"
/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/macs2 callpeak --nomodel --extsize 150 --shift -75 -t $fBam -f BAM -n $dir4"/"$RL"_NFR" --keep-dup all --gsize hs
cut -f1-3 $dir4"/"$RL"_NFR_peaks.narrowPeak" | samtools view -L - -c -b $fBam >  $dir4"/"$RL"_ReadsInNFRpeaks"
cut -f1-3 $dir4"/"$RL"_NFR_peaks.narrowPeak" | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed -u -wa -a - -b $promoters | samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInNFRpeaks_overProm"
	

# 14. MACS BAMPE peaks
echo -e "\t14. MACS2:::params for PE fragments TLEN"
/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/macs2 callpeak -t $fBam -f BAMPE -n $dir4"/"$RL"_BAMPE" --keep-dup all --gsize hs
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" | samtools view -L - -c -b $fBam >  $dir4"/"$RL"_ReadsInBAMPEpeaks"
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed -u -wa -a - -b $promoters | samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInBAMPEfrag_overProm" 
samtools view -f66 $fBam | cut -f 9 | sed 's/^-//' >  $dir4"/"$RL"_InsertSizesBAMPE"
Rscript $tool"/dependencies/Frag_hist.R" ${subBam/.bam/_InsertSizesBAMPE} $dir4"/"$RL"_InsertSizesBAMPE.png"


# 15. Reads and Peaks In Background 
echo -e "\t15. Reads in Background"
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" $dir4"/"$RL"_NFR_peaks.narrowPeak" | sort -k1,1 -k2,2n  | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/mergeBed -i - | grep -v chrM | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/bedtools subtract -wb -a $chrBed -b - | wc -l > $dir4"/"$RL"_PeaksInBackground"
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" $dir4"/"$RL"_NFR_peaks.narrowPeak" | sort -k1,1 -k2,2n  | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/mergeBed -i - | grep -v chrM | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/bedtools subtract -a $chrBed -b - | samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInBackground"

