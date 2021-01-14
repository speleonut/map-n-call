#!/bin/bash

# This is the master script that coordinates job submission for the neurogenetics WGS pipeline.
scriptDir=/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/mark/map-n-call/GATK4

usage()
{
echo "# This is the master script that coordinates job submission for primarily Illumina WGS alignments but will work for exomes too.
# The script follows GATK4x best practices as best as can be done.  
# The scripts deliver an indel realigned BAM file, a WGS metrics report, a pipeline log and a gzipped gVCF file
# Requires: BWA-MEM, samtools, sambamba, GATK, Picard, Java 
#
# Usage $0 -p file_prefix -s /path/to/sequences -o /path/to/output -c /path/to/config.cfg -S sample -L LIBRARY -I ID] | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX_R1.fastq.gz 
# -s	REQUIRED. Path to the sequence files
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory /hpcfs/users/${USER}/BWA-GATK/\$sample is used)
# -S	OPTIONAL. Sample name which will go into the BAM header and VCF header. If not specified, then will be set the same as -p
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line plain text, eg. MySeqProject20170317-PintOGuiness). Default \"IlluminaGenome\"
# -I	OPTIONAL. ID for the sequence (to go into the @RG line). If not specified the script will make one up from the first read header, and sample name
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 16/11/2017 
# Modified: (20/4/2020; Ali Gardner; for hg38 Phoenix)
# 
"
}
## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-c )			shift
					config=$1
					;;
		-p )			shift
					outPrefix=$1
					;;
		-s )			shift
					seqPath=$1
					;;
		-S )			shift
					sample=$1
					;;
		-o )			shift
					workDir=$1
					;;
		-L )			shift
					LB=$1
					;;
		-I )			shift
					ID=$1
					;;
		-h | --help )		usage
					exit 1
					;;
		* )			usage
					exit 1
	esac
	shift
done
if [ -z "$config" ]; then # If no config file specified use the default
   config=/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/mark/map-n-call/configs/BWA-GATKHC.GRCh38_full_analysis_set_phoenix.cfg
fi
if [ -z "$outPrefix" ]; then # If no file prefix specified then do not proceed
	usage
	echo "#ERROR: You need to specify a file prefix (PREFIX) referring to your sequence files eg. PREFIX_R1.fastq.gz."
	exit 1
fi
if [ -z "$seqPath" ]; then # If path to sequences not specified then do not proceed
	usage
	echo "#ERROR: You need to specify the path to your sequence files"
	exit 1
fi
if [ -z "$LB" ]; then # If library not specified then use "IlluminaGenome"
	LB=IlluminaGenome
	echo "Using $LB for library name"
fi
if [ -z "$sample" ]; then # If sample name not specified then use "outPrefix"
	sample=$outPrefix
	echo "Using $sample for sample name"
fi
if [ -z "$workDir" ]; then # If no output directory then set and create a default directory
	workDir=/hpcfs/users/${USER}/BWA-GATK/$sample
	echo "Using $workDir as the output directory"
fi
# Ensure the working directory is present or created
if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi
# Ensure the log directory is present or created
if [ ! -d "/fast/users/$USER/log" ]
    mkdir -p /fast/users/$USER/log
fi

BWAjob=`sbatch --export=ALL $scriptDir/mapSortDedupMarkIndels_hg38_phoenix.sh  -c $config -p $outPrefix -s $seqPath -S $sample -o $workDir -L $LB -I $ID`
BWAjob=$(echo $BWAjob | cut -d" " -f4)
BQSRjob=`sbatch --export=ALL --dependency=afterok:${BWAjob} $scriptDir/GATK.BQSR_hg38_phoenix.sh -c $config -p $outPrefix -S $sample -o $workDir`
BQSRjob=$(echo $BQSRjob | cut -d" " -f4)
ApplyBQSRJob=`sbatch --array=0-23 --export=ALL --dependency=afterok:${BQSRjob} $scriptDir/GATK.ApplyBQSR_hg38_phoenix.sh -c $config -p $outPrefix -S $sample -o $workDir`
ApplyBQSRJob=$(echo $ApplyBQSRJob | cut -d" " -f4)
MergeJob=`sbatch --export=ALL --dependency=afterok:${ApplyBQSRJob} $scriptDir/sambamba.Merge_hg38_phoenix.sh -c $config -p $outPrefix -S $sample -o $workDir`
MergeJob=$(echo $MergeJob | cut -d" " -f4)
GATKHCjob=`sbatch --array=0-23 --export=ALL --dependency=afterok:${MergeJob} $scriptDir/GATK.HC_hg38_phoenix.sh -c $config -p $outPrefix -S $sample -o $workDir`
GATKHCjob=$(echo $GATKHCjob | cut -d" " -f4)
sbatch --export=ALL --dependency=afterok:${MergeJob} $scriptDir/Picard.CollectWGSMetrics_hg38_phoenix.sh -c $config -p $outPrefix -S $sample -o $workDir
sbatch --export=ALL --dependency=afterok:${GATKHCjob} $scriptDir/GATK.gatherVCFs_hg38_phoenix.sh -c $config -p $outPrefix -S $sample -o $workDir
