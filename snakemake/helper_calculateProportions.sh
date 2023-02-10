#!/bin/bash


numberOfReads="$1" # Specificed by USER: arg.1
basename="$2" # Specificed by USER: arg.2
output="$3" # Specificed by USER: arg.3


reads=$(cat $numberOfReads) 

printf '%s\n' {1000000..10000000..1000000}/$reads | bc -l | xargs printf '%.3f\n' | awk -v awkbase=$basename '$1<1 {print awkbase " " $0}' > $output
printf '%s\n' {15000000..40000000..5000000}/$reads | bc -l | xargs printf '%.3f\n' | awk -v awkbase=$basename '$1<1 {print awkbase " " $0}' >> $output
