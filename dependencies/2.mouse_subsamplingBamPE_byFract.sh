#!/bin/bash
######################################################################
# Subsamplig a bam file base on fractions previously calculated
#
# Creates several output files with all the metrics: 
#
# Example: 
#./3.subsamplingBamPE_byFract.sh /scratchfs/aalvarez/Oliver_BulkATAC/Outputs_Novaseq/4.Alignment/RL1753.bam > RL1753_SubsSamp.log 2>&1 &
#
# This script is included in the subsampling pipeline !!!
######################################################################

# If any command fails or a pipe breaks, the script will stop running.
set -e -o pipefail -o verbose

##### Inputs:

bam="$1" # Specificed by USER: arg.1
dirO="$2" # Specificed by USER: arg.2
tool="/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool" # Specificed by USER: arg.3

RL=$(basename ${bam/.bam/})

# Some utilis:
sambamba=/home/sbuckberry/working_data_01/bin/sambamba_v0.6.3
parallel=/home/sbuckberry/working_data_01/bin/parallel
blacklist=$tool"/data/ENCFF547MET.bed"
chromSizes=$tool"/data/mm10.chrom.sizes"
promoters="/home/dmakosa/working_data_04/genomes/mouse/ucsc/mm10/epd_newpromoter/epd_newpromoter_upstream2000.bed"
chrBed=$tool"/data/mm10.chromS.bed"
black_M=$tool"/data/mm10.chrM.size.bed"
allRLsNFR=$tool"/data/RL2202_20M_pairs_dedup_filt_noMT_NFR_peaks.narrowPeak"
Rhist=$tool"/dependencies/Frag_hist.R"
MACS2="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/macs2"
intBed="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed"
mergeBed="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/mergeBed"
BedTools="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/bedtools"

# Output:
dirOut="$dirO/5.Subsampling/$RL"

echo -e "\t\tWorking with bam: $bam in $dirOut\n\t\tID: $RL\n\n"  

# catch proportions:
pro=$dirOut"/"$RL".proportions"
se="0 + 0 read2" 
seed=39 # seed for sambamba... reproducibility

cut -f2 $pro | tail -n+3 | while read p # loop through proportions of interest
do
	echo -e "\t1. Subsampling proportion of "$p
	subBam=$dirOut"/"$RL"_"$p".bam" #subOut
	  # subsampling
	$sambamba view -h -t 4 -s $p -f bam --subsampling-seed $seed $bam -o $subBam
	
	  # stats check subsam ok
	samtools flagstat $subBam > ${subBam/.bam/_flagStats}	
	
	  # 2) Auto check PE or SE, proper pairs
	pairs=$(grep read2 ${subBam/.bam/_flagStats})
	if [ "$pairs" == "$se" ] #SE
	then 
		echo -e "\t2. Data is SE --> $pairs"
		mv $subBam ${subBam/.bam/merged_pairs.bam}
	else
		echo -e "\t2. Data is PE ---> $pairs\n\t\tFilter proper Pairs" 
		$sambamba view -t 4 -f bam -F "proper_pair" $subBam -o ${subBam/.bam/merged_pairs.bam}
		samtools flagstat ${subBam/.bam/merged_pairs.bam} > ${subBam/.bam/merged_pairs_flagStats}
	fi

	   # 3 Count MT reads
	echo -e "\t3. Count MT reads"
	samtools index ${subBam/.bam/merged_pairs.bam} # check if needed!
	samtools view -c ${subBam/.bam/merged_pairs.bam} chrM > ${subBam/.bam/merged_pairs_mitoReads}	

	
	   # 4 PCR dups
	echo -e "\t4. Remove PCR dups"  
	$sambamba markdup -p -t 4 -r --hash-table-size=1000000 ${subBam/.bam/merged_pairs.bam} ${subBam/.bam/merged_pairs_dedup.bam}
	samtools flagstat ${subBam/.bam/merged_pairs_dedup.bam} > ${subBam/.bam/merged_pairs_dedup_flagStats}
	rm ${subBam/.bam/merged_pairs.bam}	
	rm ${subBam/.bam/merged_pairs.bam.bai}

	  # 5 BlackList Regions + Mito Reads 
	echo -e "\t5. Intercept BlackListed regions and Mito reads"
	samtools view -q 30 -L $blacklist -U ${subBam/.bam/merged_pairs_dedup_filt.bam} -b ${subBam/.bam/merged_pairs_dedup.bam} > ${subBam/.bam/_crap}
	rm ${subBam/.bam/merged_pairs_dedup.bam}

	samtools view -L $black_M -U ${subBam/.bam/merged_pairs_dedup_filt_noMT.bam} -b ${subBam/.bam/merged_pairs_dedup_filt.bam}  >  ${subBam/.bam/_2crap} 

	samtools flagstat ${subBam/.bam/merged_pairs_dedup_filt.bam} > ${subBam/.bam/merged_pairs_dedup_filt_flagStats} 
	rm ${subBam/.bam/merged_pairs_dedup_filt.bam} 
	samtools flagstat ${subBam/.bam/merged_pairs_dedup_filt_noMT.bam} > ${subBam/.bam/merged_pairs_dedup_filt_noMT_flagStats} 
	
	  # 6 Rm intermediate
	echo -e "\t6. Remove Intermediate Files"
	rm $subBam
	rm $subBam".bai"
	rm ${subBam/.bam/merged_pairs_dedup.bam.bai}
	rm ${subBam/.bam/}*crap
	
	  # redefine var
	fBam=${subBam/.bam/merged_pairs_dedup_filt_noMT.bam}

	  # 7 Reads in Promoters:
	echo -e "\t 7. Calcule Reads in Promoters"
	# samtools counts reads or pairs? 
	# samtools view -L $promoters  $fBam | wc -l  --> out is same as Reads in Prom , therefore for PE is counting: PAIRs
	samtools view -c -L $promoters -b $fBam > ${subBam/.bam/_ReadsInProm}

	  #  8 and 9 MACS2 : Paired-End mode is off just check MACS2 out
	echo -e "\t 8. MACS2:::params for Tn5 insert sites, params ENCODE pipe for NFRs"
	$MACS2 callpeak --nomodel --extsize 150 --shift -75 -t $fBam -f BAM -n ${subBam/.bam/_NFR} --keep-dup all --gsize mm	

	#echo -e "\t 9. MACS2:::params for Tn5 insert sites, 200 and 100 "
        # NOPE: macs2 callpeak --nomodel --extsize 200 --shift -100 -t $fBam -f BAM -n ${subBam/.bam/_Tn5Insert} --keep-dup all --gsize hs

	# Reads in Peaks: what are the peaks I want? fixed size or use narrowOutput? (narrow)
	echo -e "\t 10. Calcule Reads in NFR Peaks"
	cut -f1-3 ${subBam/.bam/_NFR_peaks.narrowPeak} | samtools view -L - -c -b $fBam >  ${subBam/.bam/_ReadsInNFRpeaks}
	# NOPE: cut -f1-3 ${subBam/.bam/_Tn5Insert_peaks.narrowPeak} | samtools view -L - -c -b $fBam > ${subBam/.bam/_ReadsInTn5peaks}
	
	# Reads Peaks overlapping promotors
	echo -e "\t 11. Calcule Reads in NFR Peaks overlapping a promoter"
	cut -f1-3 ${subBam/.bam/_NFR_peaks.narrowPeak} | $intBed -u -wa -a - -b $promoters | samtools view -L - -c -b $fBam > ${subBam/.bam/_ReadsInNFRpeaks_overProm}
	# NOPE: cut -f1-3 ${subBam/.bam/_Tn5Insert_peaks.narrowPeak} | intersectBed -wa -a - -b $promoters | samtools view -L - -c -b $fBam >  ${subBam/.bam/_ReadsInTn5peaks_overProm}

	# Calcule Reads in NFR peaks of the track with All RLs from Novaseq experiment merged:
		# experiments here: /dd_rundata/novaseq/Runs/190608_A00690_H7V22DRXX_004/Unaligned/snATAC/RL*
	echo -e "\t 12. Calcule Reads in Peaks all RLs "
        cut -f1-3 $allRLsNFR | samtools view -L - -c -b $fBam >  ${subBam/.bam/_ReadsInNFRpeaks_fromALL_RLsNovaseq}

	#  only 4 PE
	if [ "$pairs" == "$se" ] #SE
	then 
		echo -e "\tData is SE --> NO BAMPE"
		# Reads in background: join all 2 sets of peaks and calculate reads out of peaks:
		echo -e "\t 13. Reads in Background"
		cut -f1-3 ${subBam/.bam/_NFR_peaks.narrowPeak} ${subBam/.bam/_Tn5Insert_peaks.narrowPeak} | sort -k1,1 -k2,2n  | $mergeBed -i - | grep -v chrM | $BedTools subtract -a $chrBed -b - | samtools view -L - -c -b $fBam > ${subBam/.bam/_ReadsInBackground}
	else
		echo -e "\tData is PE ---> BAMPE option" 
		echo -e "\t13. MACS2:::params for PE fragments TLEN"
		$MACS2 callpeak -t $fBam -f BAMPE -n ${subBam/.bam/_BAMPE} --keep-dup all --gsize mm	
		echo -e "\t14 Reads in Peaks"
		cut -f1-3 ${subBam/.bam/_BAMPE_peaks.narrowPeak} | samtools view -L - -c -b $fBam >  ${subBam/.bam/_ReadsInBAMPEfrag}
		echo -e "\t 15. Reads in Peaks over promotors"
		cut -f1-3 ${subBam/.bam/_BAMPE_peaks.narrowPeak} | $intBed -u -wa -a - -b $promoters | samtools view -L - -c -b $fBam > ${subBam/.bam/_ReadsInBAMPEfrag_overProm}	
		# Reproduce histogram from Bioanalizer: histogram of insertSizes: use phantom and also repro what MACS is doing with no Model
		# 66 flag: mapped in proper pair and first in pair: only PE guy ; TLEN is field 9 (some are negatives)
		echo -e "\t 16. Insert sizes BAMPE and Plot histogram"
		samtools view -f66 $fBam | cut -f 9 | sed 's/^-//' >  ${subBam/.bam/_InsertSizesBAMPE}
		Rscript $Rhist ${subBam/.bam/_InsertSizesBAMPE} ${subBam/.bam/_InsertSizesBAMPE.png}

		echo -e "\t 17. Reads in Background"
		# NOPE: cut -f1-3 ${subBam/.bam/_NFR_peaks.narrowPeak} ${subBam/.bam/_Tn5Insert_peaks.narrowPeak} ${subBam/.bam/_BAMPE_peaks.narrowPeak} | sort -k1,1 -k2,2n  | mergeBed -i - | grep -v chrM | bedtools subtract -a $chrBed -b - | samtools view -L - -c -b $fBam > ${subBam/.bam/_ReadsInBackground}
		cut -f1-3 ${subBam/.bam/_NFR_peaks.narrowPeak} ${subBam/.bam/_BAMPE_peaks.narrowPeak} | sort -k1,1 -k2,2n  | $mergeBed -i - | grep -v chrM | $BedTools subtract -a $chrBed -b - | samtools view -L - -c -b $fBam > ${subBam/.bam/_ReadsInBackground}
	fi
	
	# rm bam file 	
	echo -e "\t\tbye bye bam\n\n\n"
	rm $fBam

	#break
done



