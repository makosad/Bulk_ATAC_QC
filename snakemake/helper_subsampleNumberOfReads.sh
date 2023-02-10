#!/bin/bash


proportions="$1" # Specificed by USER: arg.1
basename="$2" # Specificed by USER: arg.2
output="$3" # Specificed by USER: arg.3

export PROMOTERS=$promoters
cat $proportions | xargs -n2 sh -c 'samtools view -c output/5.SubsamplingOfFiltered/1.subsampling/sub$1_$0_pairs_dedup_filt_noMT.bam > output/5.SubsamplingOfFiltered/3.readsInPeaks/sub$1_$0_pairs_dedup_filt_noMT.NoOfReads'
touch $output

