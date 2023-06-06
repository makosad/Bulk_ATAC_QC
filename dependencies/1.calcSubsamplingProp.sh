#!/bin/bash
###################################################################################
# Defining subsampling fractions based on TOTAL READs (_flagSTat files created previously)
#
# Creates output txt file with numbers: Total Reads and proportion for subsampling
#
# To run this pipe with your own samples, please just provide the corresponding 2 arguments: 
# ("$1") path to read1 and ("$2") path to output folder
# they have tag at the end:     # Specificed by USER
#
# Example run: 
#	./2.calcSubsamplingProp.sh /home/aalvarez/Bulk_ATAC_QCmetrics_tool/output/4.Alignment/RL1753_pairs_dedup_filt_noMT.bam /home/aalvarez/Bulk_ATAC_QCmetrics_tool/output/
#
#
# This script is included in the subsampling pipeline !!!  
###################################################################################

# If any command fails or a pipe breaks, the script will stop running.
set -e -o pipefail -o verbose

###### Inputs:

bam="$1" # Specificed by USER: arg.1

dO="$2" # Specificed by USER: arg.2
dOut=${dO}"/5.Subsampling"

mkdir -p $dOut

RL=$(basename ${bam/.bam/})
TR=$(grep QC ${bam/.bam/_flagStats} | sed 's/ .*//g') 

mkdir -p $dOut"/"$RL

	# open FILE
echo -e "#TotalReads=$TR\nReads\tProportion"  >  $dOut/$RL/$RL.proportions	

# for f in {1..10} 20 25 30 35 40 45 50 60 70 80 90 100 110 120 130 140 150 # calcule proportion for these numbers of Million reads !! Change if desired

for f in {1..10} 20 25 30 35 100  
do 
	f=$(($f*1000000)) #million reads
	pro=$(echo "scale=3 ; $f/$TR"| bc ) # desired proportion
	# if pro > 1 avoid
	echo -e $f"\t"$pro | awk '$2 < 1 {print $0}'  >> $dOut/${RL}/${RL}.proportions
done

