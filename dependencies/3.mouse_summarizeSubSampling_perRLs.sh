#!/bin/bash

# If any command fails or a pipe breaks, the script will stop running.
set -e -o pipefail -o verbose

bam="$1"
dirO="$2"

RL=$(basename ${bam/.bam/})
dirI=$dirO"/5.Subsampling/"$RL
dirS=${dirI}"/StatsAndPlots"
promoters="/home/dmakosa/working_data_04/genomes/mouse/ucsc/mm10/epd_newpromoter/epd_newpromoter_upstream2000.bed"
intBed="/home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed"
mkdir -p $dirS

echo "#" $(date) > ${dirS}/${RL}_SubSamplingSummary.txt # start file
echo -e "RL\tproportion\tTotalReads\tMapped\tFilt\tReadsInProm\tNFRpeaks\tReadsInNFRpeaks\tReadsInNFRpeaksOverProm\tNFRPeaksInPromoters\tReadsInBackground\tBAMPEPeaks\tReadsInBAMPEpeaks\tReadsInBAMPEpeaksOverProm\tBAMPEPeaksInPromoters\treadsin_NFRPeaksAllRLs" >> ${dirS}/${RL}_SubSamplingSummary.txt

for p in $(ls -rt $dirI/*NFR_peaks.xls) # timesort of fractions
do
	p=$(basename ${p/_NFR*/} | sed "s/.*_//")
	echo $p

	#TOTAL"sequenced"reads (no Pairs)
	uno=$(grep QC $(ls $dirI/*${p}_flagStats ) | sed 's/+.*//g' )

	dos=$(grep "mapped ("  $(ls $dirI/*${p}_flagStats ) | sed 's/+.*//g' )
	#onlyPE:dos=$(grepread1$(ls-tr$name"_"$f*Stats|tail-n3)|sed's/+.*//g'|sed's/.*://g'|tr'\n''\t')
		
	validR=$(grep QC  $dirI/*${p}*noMT_flagStats |  sed 's/+.*//g')

	#Reads(orFragments!!!)inPromotors
	tres=$(cat $dirI/*_$p*_ReadsInProm )

	  # NFR peaks
	#TotalnumberofNFRPeaks
	cuatro=$(wc -l $dirI/*_$p*NFR_summits.bed | sed 's/ .*//')
	#ReadsinNFRPeaks
	cinco=$(cat $dirI/*_$p*_ReadsInNFRpeaks)	
	#ReadsinNFRPeaksOverPromotors
	seis=$(cat $dirI/*_$p*_ReadsInNFRpeaks_overProm)
	#NumberPeaksinNFRoverlappingpromotors
	siete=$(cut -f1-3 $dirI/*_$p*NFR_peaks.narrowPeak | $intBed -u -wa -a - -b $promoters | wc -l)	
		
	#TotalnumberofTN5Peaks
	#ocho=$(wc -l $dirI/*_$p*_Tn5Insert_summits.bed | sed 's/ .*//')
	#ReadsinTN5Peaks
	#nueve=$(cat $dirI/*_$p*_ReadsInTn5peaks)
	#ReadsinTN5PeaksOverPromotors
	#diez=$(cat $dirI/*_$p*_ReadsInTn5peaks_overProm)
	#NumberPeaksinNFRoverlappingpromotors#idiotshoulddidthisinscript3
	#once=$(cut -f1-3 $dirI/*_$p*_Tn5Insert_peaks.narrowPeak | intersectBed -wa -a - -b /scratchfs/aalvarez/Oliver_BulkATAC/promoterhg19_2kb.bed | wc -l)	

	#ReadsinBackground		
	doce=$(cat $dirI/*_$p*_ReadsInBackground)
	
 	  # BAMPE-->onlyPE
	#TotalnumberofBAMPEPeaks
	trece=$(wc -l $dirI/*_$p*_BAMPE_summits.bed | sed 's/ .*//')
	#ReadsinBAMPEfrags
	catorce=$(cat $dirI/*_$p*_ReadsInBAMPEfrag)
	#ReadsinBAMPEOverPromotors
	quince=$(cat $dirI/*_$p*_ReadsInBAMPEfrag_overProm)
	#NumberPeaksinBAMPEoverlappingpromotors#idiotshoulddidthisinscript3
	dieciseis=$(cut -f1-3 $dirI/*_$p*BAMPE_peaks.narrowPeak | $intBed -u -wa -a - -b $promoters | wc -l)	
	#OUTFILEperRL

	allRLs=$(cat $dirI/*_$p*_ReadsInNFRpeaks_fromALL_RLsNovaseq)

	echo -e "$RL\t$p\t$uno\t$dos\t$validR\t$tres\t$cuatro\t$cinco\t$seis\t$siete\t$doce\t$trece\t$catorce\t$quince\t$dieciseis\t$allRLs" >>  ${dirS}/${RL}_SubSamplingSummary.txt
	#>>$dirOut"/all_RL_stats_v2.tsv"
	# sed-i's///g' allStats_v2.txt

	#break
done


