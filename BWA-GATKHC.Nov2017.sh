#!/bin/bash

# This is the master script that coordinates job submission for the nuerogenetics WGS pipeline.
scriptDir=/data/neurogenetics/git/PhoenixScripts/shared/scripts/

usage()
{
echo "# This is the master script that coordinates job submission for primarily Illumina WGS alignments but will work for exomes too.
# The script is considerably restructured and hopefully will now suit the architecture of Phoenix.  
# The scripts deliver an indel realigned BAM file, a WGS metrics report, a pipeline log and a gzipped gVCF file
# Requires: BWA-MEM, samtools, sambamba, GATK, Picard, Java 
#
# Usage $0 -p file_prefix -s /path/to/sequences -o /path/to/output -c /path/to/config.cfg -S SAMPLE -L LIBRARY -I ID] | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX_R1.fastq.gz 
# -s	REQUIRED. Path to the sequence files
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $FASTDIR/BWA-GATK/\$sample is used)
# -S	OPTIONAL. Sample name which will go into the BAM header and VCF header. If not specified, then will be set the same as -p
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line plain text, eg. MySeqProject20170317-PintOGuiness). Default \"IlluminaGenome\"
# -I	OPTIONAL. ID for the sequence (to go into the @RG line). If not specified the script will make one up from the first read header, and sample name
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 16/11/2017 
# Modified: (Date; Name; Description)
# 16/04/2020; Mark Corbett; Fix quoting, fix log directory, remove redundant information, fix variable names
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
   config=/data/neurogenetics/git/PhoenixScripts/shared/scripts/BWA-GATKHC.genome.cfg
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
	echo "Using $outPrefix for sample name"
fi
if [ -z "$workDir" ]; then # If no output directory then set and create a default directory
	workDir=$FASTDIR/BWA-GATK/$sample
	echo "Using $FASTDIR/BWA-GATK/$sample as the output directory"
fi
if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi
# Ensure the log directory is present or created
if [ ! -d "/fast/users/$USER/log" ]
    mkdir -p /fast/users/$USER/log
fi

# Submit jobs to the queue
BWAjob=`sbatch --export=ALL $scriptDir/mapSortDedupMarkIndels.sh  -c $config -p $outPrefix -s $seqPath -S $sample -o $workDir -L $LB -I $ID`
BWAjob=$(echo $BWAjob | cut -d" " -f4)
InDelJob=`sbatch --array=0-23 --export=ALL --dependency=afterok:${BWAjob} $scriptDir/GATK.IndelRealign.sh -c $config -p $outPrefix -S $sample -o $workDir`
InDelJob=$(echo $InDelJob | cut -d" " -f4)
BQSRjob=`sbatch --export=ALL --dependency=afterok:${InDelJob} $scriptDir/GATK.BQSR.sh -c $config -p $outPrefix -S $sample -o $workDir`
BQSRjob=$(echo $BQSRjob | cut -d" " -f4)
PrintReadsJob=`sbatch --array=0-23 --export=ALL --dependency=afterok:${BQSRjob} $scriptDir/GATK.PrintReads.sh -c $config -p $outPrefix -S $sample -o $workDir`
PrintReadsJob=$(echo $PrintReadsJob | cut -d" " -f4)
MergeJob=`sbatch --export=ALL --dependency=afterok:${PrintReadsJob} $scriptDir/sambamba.Merge.sh -c $config -p $outPrefix -S $sample -o $workDir`
MergeJob=$(echo $MergeJob | cut -d" " -f4)
GATKHCjob=`sbatch --array=0-23 --export=ALL --dependency=afterok:${MergeJob} $scriptDir/GATK.HC.sh -c $config -p $outPrefix -S $sample -o $workDir`
GATKHCjob=$(echo $GATKHCjob | cut -d" " -f4)
sbatch --export=ALL --dependency=afterok:${MergeJob} $scriptDir/Picard.CollectWGSMetrics.sh -c $config -p $outPrefix -S $sample -o $workDir
sbatch --export=ALL --dependency=afterok:${GATKHCjob} $scriptDir/GATK.VCFconcat.sh -c $config -p $outPrefix -S $sample -o $workDir
