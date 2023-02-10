#!/bin/bash


proportions="$1" # Specificed by USER: arg.1
basename="$2" # Specificed by USER: arg.2
output="$3" # Specificed by USER: arg.3
promoters="$4" # Specificed by USER: arg.4

export PROMOTERS=$promoters
cat $proportions | xargs -n2 sh -c 'cut -f1-3 output/5.SubsamplingOfFiltered/2.callpeaks/sub$1_$0_pairs_dedup_filt_noMT.NFR_peaks.narrowPeak | samtools view -L - -c -b output/5.SubsamplingOfFiltered/1.subsampling/sub$1_$0_pairs_dedup_filt_noMT.bam > output/5.SubsamplingOfFiltered/3.readsInPeaks/sub$1_$0_pairs_dedup_filt_noMT.ReadsInNFRpeaks'
cat $proportions | xargs -n2 sh -c 'cut -f1-3 output/5.SubsamplingOfFiltered/2.callpeaks/sub$1_$0_pairs_dedup_filt_noMT.NFR_peaks.narrowPeak | /home/dmakosa/working_data_01/apps/miniconda3/envs/deeptoolsenv/bin/intersectBed -u -wa -a - -b "$PROMOTERS" | samtools view -L - -c -b output/5.SubsamplingOfFiltered/1.subsampling/sub$1_$0_pairs_dedup_filt_noMT.bam > output/5.SubsamplingOfFiltered/3.readsInPeaks/sub$1_$0_pairs_dedup_filt_noMT.ReadsInNFRpeaks_overProm'
touch $output


