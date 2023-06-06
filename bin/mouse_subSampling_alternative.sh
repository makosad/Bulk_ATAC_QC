#!/bin/bash


######################################################################
#####    READMEEEEEE !!!! 
######################################################################
# 	Subsamplig a bam file and perform bulkATACseq analysis:
#
#
	# Edit this file ONLY !!
#
# IMPORTANT: 
 
	# _flagStat file must be in the same dir than bam file !!!
	#  if not, create!!
#
# samtools flagstat $bam > ${bam/.bam/_flagStats} # _flagStat file with metrics
#
#
# ADVICEs and TIPs: 
#
#	- use the SAME output dir ($dirO arg.2) you used with the analyzeBulkATACseqRLs_AllReads.sh step
#	- if you are gonne to subsample several bamfiles, use the SAME output dir ($dirO arg.2) for all of them 
#	- maximum CPUs reached by this script --> 5 (just at some steps)
#	- if cp/paste this Tool to your home or whatever, cp the FULL thing:
#		  cp -r Bulk_ATAC_QCmetrics_tool/ yourFavDir/.
#	  and then change the correspondin path --> tool="" 
#	
######################################################################

# If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

####### Inputs to Specify by the USER:

# alignment file: the original one, NO the processed one!!!
bam="$1"     # Specificed by USER: FULL PATH (exmpl: /scratchfs/aalvarez/Oliver_BulkATAC/Outputs_Novaseq/4.Alignment/RL1753.bam) 

# output main dir
dirO="$2"    # Specificed by USER: FULL PATH (exmpl: /scratchfs/aalvarez/Oliver_BulkATAC/Outputs_Novaseq )

tool="/home/dmakosa/working_data_04/Bulk_ATAC_QCmetrics_tool"  # Specificed by USER: FULL PATH (exmpl: /home/aalvarez/Bulk_ATAC_QCmetrics_tool)


####### Run pipeline: 

echo -e "\t\tStep 1:: cacule subsampling fractions\n"
#${tool}/dependencies/1.calcSubsamplingProp.sh $bam $dirO  # works fine


echo -e "\t\tStep 2:: subsampling bam file \n"
#${tool}/dependencies/2.mouse_subsamplingBamPE_byFract.sh $bam $dirO $tool # works

echo -e "\t\tStep 3:: summarize subsampling data per RL\n"
#${tool}/dependencies/3.mouse_summarizeSubSampling_perRLs.sh $bam $dirO  # works

echo -e "\t\tStep 4:: plots per RL\n"
RL=$(basename ${bam/.bam/})
txt=$dirO"/5.Subsampling/${RL}/StatsAndPlots/${RL}_SubSamplingSummary.txt"
Rscript ${tool}/dependencies/4.alternative.R $txt 

# echo -e "\t\tStep 5:: \n"

