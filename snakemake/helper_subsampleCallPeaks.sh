#!/bin/bash


proportions="$1" # Specificed by USER: arg.1
basename="$2" # Specificed by USER: arg.2
output="$3" # Specificed by USER: arg.3
size="$4" # Specificed by USER: arg.4


export SIZE=$size
cat $proportions | xargs -n2 sh -c 'macs2 callpeak --nomodel --extsize 150 --shift -75 -t output/5.SubsamplingOfFiltered/1.subsampling/sub$1_$0_pairs_dedup_filt_noMT.bam -f BAM -n output/5.SubsamplingOfFiltered/2.callpeaks/sub$1_$0_pairs_dedup_filt_noMT.NFR --keep-dup all --gsize "$SIZE" 2> output/5.SubsamplingOfFiltered/2.callpeaks/log_sub$1_$0_pairs_dedup_filt_noMT.log'
touch $output
