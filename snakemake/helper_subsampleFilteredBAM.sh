#!/bin/bash


proportions="$1" # Specificed by USER: arg.1
basename="$2" # Specificed by USER: arg.2
output="$3" # Specificed by USER: arg.3
threads="$4" # Specificed by USER: arg.4

export THREADS=$threads
cat $proportions | xargs -n2 sh -c '/home/sbuckberry/working_data_01/bin/sambamba_v0.6.3 view -h -t "$THREADS" -s $1 -f bam --subsampling-seed 39 output/4.Alignment/mapped_pairs_dedup_filt_noMT/$0_pairs_dedup_filt_noMT.bam -o output/5.SubsamplingOfFiltered/1.subsampling/sub$1_$0_pairs_dedup_filt_noMT.bam'
touch $output

