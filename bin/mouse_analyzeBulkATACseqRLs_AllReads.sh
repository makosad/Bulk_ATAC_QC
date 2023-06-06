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
# sh ~/working_data_04/Bulk_ATAC_QCmetrics_tool/bin/mouse_analyzeBulkATACseqRLs_AllReads.sh /scratchfs/dmakosa/bulk10xmethylomes/210212/data/RL2384new_2021_02_12_methanolNuclei_lambda_ezdna_S19_R1_001.fastq.gz /scratchfs/dmakosa/bulk10xmethylomes/210212/analyzeasbulkatac > RL2384new.log  2>&1 &  
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
set -e -o pipefail -o verbose

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
FASTQC="/home/dmakosa/working_data_01/apps/miniconda3/envs/10xmethylomes/bin/fastqc"
#BBDUK="/home/sbuckberry/working_data_01/bin/bbmap/bbduk2.sh"
CUTADAPT="/home/dmakosa/working_data_01/apps/miniconda3/envs/10xmethylomes/bin/cutadapt"
SAMBAMBA="/home/sbuckberry/working_data_01/bin/sambamba_v0.6.3"
parallel="/home/dmakosa/working_data_01/apps/miniconda3/envs/10xmethylomes/bin/parallel"
blacklist="/home/dmakosa/working_data_04/genomes/mouse/ucsc/mm10/Blacklistedregions/ENCFF547MET.bed"
chromSizes="/home/dmakosa/working_data_04/genomes/mouse/ucsc/mm10/mm10.chrom.sizes"
promoters="/home/dmakosa/working_data_04/genomes/mouse/ucsc/mm10/Downloaded_additions/Promoters_UCSC_GenesAndPredictions_allGENCODE_VM25_Basic.bed" 
chrBed=$tool"/data/mm10.chromS.bed"
index="/home/dmakosa/working_data_04/genomes/mouse/ucsc/mm10/UCSC_mm10_genome_Bowtie2Index/mm10_genome"
black_M="/home/dmakosa/working_data_04/genomes/mouse/ucsc/mm10/mm10.chrM.size.bed"
MACS2="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/macs2"
intBed="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed"
mergeBed="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/mergeBed"
BedTools="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/bedtools"
bamCOV="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/bamCoverage"


###### Aanalyze:

# 1.QCs raw reads
echo -e "\t1. QCs"
$FASTQC -o $dir1 $R1 $R2 > ${dir0}/logPreQCs_${RL}

# 2.Clean Adaptors	
echo -e "\t2. Clean Adaptors"
$CUTADAPT -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o $dir2"/"$RL"_R1_clean.fq.gz" -p $dir2"/"$RL"_R2_clean.fq.gz" $R1 $R2 > ${dir0}"/logAdap"$RL
	
# 3. QCs after clean reads
echo -e "\t3. QCs after adaptors"
$FASTQC -o $dir3  $dir2"/"$RL"_R1_clean.fq.gz"  $dir2"/"$RL"_R2_clean.fq.gz" > ${dir0}/logQCs_${RL}

#4. Bowtie2 alignment::care or -T option, do not overlap temp files
echo -e "\t4. Bowtie and samtools sort" 
(bowtie2 -q --threads 8 -X2000 -x $index -1 $dir2"/"$RL"_R1_clean.fq.gz"  -2 $dir2"/"$RL"_R2_clean.fq.gz" | samtools view -bSu - | samtools sort -T ${RL}sorted - > $dir4"/"${RL}.bam) >  ${dir0}/logBowtie${RL}

# 5. STATS
echo -e "\t5. Flag Stats"
samtools flagstat $dir4"/"${RL}.bam  > $dir4"/"${RL}_flagStats

# 6. Count MT reads
echo -e "\t6. Count MT reads"
samtools index  $dir4"/"${RL}.bam # needed for -c
samtools view -c $dir4"/"${RL}.bam chrM > $dir4"/"${RL}_mitoReads


# 7. Proper Pairs
echo -e "\t7. Data is PE ---> $pairs\n\tFilter alignment, proper Pairs and stats" 
$SAMBAMBA view -t 4 -f bam -F "proper_pair" $dir4"/"${RL}.bam -o $dir4"/"${RL}_pairs.bam 
samtools flagstat $dir4"/"${RL}_pairs.bam > $dir4"/"${RL}_pairs_flagStats


# 8. PCR dups
echo -e "\t8. Remove PCR dups and stats"
$SAMBAMBA markdup -p -t 4 -r --hash-table-size=1000000  $dir4"/"${RL}_pairs.bam $dir4"/"${RL}_pairs_dedup.bam
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
$bamCOV --binSize 1 -p 5 -b $dir4"/"${RL}_pairs_dedup_filt_noMT.bam -o $dir4"/"${RL}_pairs_dedup_filt_noMT_coverage.bigwig # maybe problems here, depending on size

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
$MACS2 callpeak --nomodel --extsize 150 --shift -75 -t $fBam -f BAM -n $dir4"/"$RL"_NFR" --keep-dup all --gsize mm
cut -f1-3 $dir4"/"$RL"_NFR_peaks.narrowPeak" | samtools view -L - -c -b $fBam >  $dir4"/"$RL"_ReadsInNFRpeaks"
cut -f1-3 $dir4"/"$RL"_NFR_peaks.narrowPeak" | $intBed -u -wa -a - -b $promoters | samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInNFRpeaks_overProm"
	

# 14. MACS BAMPE peaks
echo -e "\t14. MACS2:::params for PE fragments TLEN"
$MACS2 callpeak -t $fBam -f BAMPE -n $dir4"/"$RL"_BAMPE" --keep-dup all --gsize mm
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" | samtools view -L - -c -b $fBam >  $dir4"/"$RL"_ReadsInBAMPEpeaks"
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" | $intBed -u -wa -a - -b $promoters | samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInBAMPEfrag_overProm" 
samtools view -f66 $fBam | cut -f 9 | sed 's/^-//' >  $dir4"/"$RL"_InsertSizesBAMPE"
#Rscript /home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool/dependencies/Frag_hist.R ${subBam/.bam/_InsertSizesBAMPE} $dir4"/"$RL"_InsertSizesBAMPE.png"
#Commented, because it's not working and breaks the pipeline... will investigate later

# 15. Reads and Peaks In Background 
echo -e "\t15. Reads in Background"
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" $dir4"/"$RL"_NFR_peaks.narrowPeak" | sort -k1,1 -k2,2n  | $mergeBed -i - | grep -v chrM | $BedTools subtract -wb -a $chrBed -b - | wc -l > $dir4"/"$RL"_PeaksInBackground"
cut -f1-3 $dir4"/"$RL"_BAMPE_peaks.narrowPeak" $dir4"/"$RL"_NFR_peaks.narrowPeak" | sort -k1,1 -k2,2n  | $mergeBed -i - | grep -v chrM | $BedTools subtract -a $chrBed -b - | samtools view -L - -c -b $fBam > $dir4"/"$RL"_ReadsInBackground"



