#!/bin/bash

#SBATCH -J ApplyBQSR
#SBATCH -o /hpcfs/users/%u/launch/slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=05:00:00
#SBATCH --mem=8GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# load modules
module load arch/haswell
module load Java/1.8.0_121
### module load GATK 4.0.0.0
### module load picard/2.6.0 or higher

# run the executable
# A script to map reads and then call variants using the GATK v4.x best practices designed for the Phoenix supercomputer but will work on stand alone machines too

# Variables that usually don't need changing once set for your system

usage()
{
echo "# Script for processing and mapping Illumina 100bp pair-end sequence data and optionally plotting coverage for an interval
# Requires: BWA 0.7.x, Picard, samtools, GATKv4.x, BWA-Picard-GATK-CleanUp.sh.
# This script assumes your sequence files are gzipped
#
# Usage $0 -p file_prefix -s /path/to/sequences -o /path/to/output [-i /path/to/bedfile.bed] | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX_R1.fastq.gz 
# -S	REQUIRED. Sample name if not specified then will be set the same as -p
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified current directory is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# 
# Original: Derived from Illumina-Phred33-PE-FASTX-BWA-Picard-GATKv2.sh by Mark Corbett, 17/03/2014
# Modified: (Date; Name; Description)
# 24/09/2015; Mark Corbett; Fork original for genomes
# 25/09/2015; Mark Corbett; Pipe to samtools sort; Add getWGSMetrics
# 12/10/2015; Mark Corbett; Fix error collecting .bam files to merge
# 13/05/2016; Mark Corbett; Add option to specify sample name different from OUTPREFIX.  Make seq file search explicit for *.fastq.gz
# 01/07/2016; Mark Corbett; Improve error handling
# 24/08/2016; Mark Corbett; Fork for HPC version, bring up to date with GATKv3.6
# 18/11/2016; Mark Corbett; Step down number of splits for PrintReads for higher efficiency
# 17/05/2017; Atma Ivancevic; translating for SLURM
# 16/11/2017; Mark Corbett; Fork to split scripts to see if this works on Phoenix; swap Picard for sambamba for duplicate marking
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-c )			shift
					CONFIG=$1
					;;
		-p )			shift
					OUTPREFIX=$1
					;;
		-S )			shift
					SAMPLE=$1
					;;
		-o )			shift
					WORKDIR=$1
					;;
		-h | --help )		usage
					exit 0
					;;
		* )			usage
					exit 1
	esac
	shift
done
if [ -z "$CONFIG" ]; then # If no config file specified use the default
   CONFIG=/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/GATK4/BWA-GATKHC.GRCh37_full_analysis_set_phoenix.cfg
fi
source $CONFIG

if [ -z "$OUTPREFIX" ]; then # If no file prefix specified then do not proceed
	usage
	echo "#ERROR: You need to specify a file prefix (PREFIX) referring to your sequence files eg. PREFIX_R1.fastq.gz."
	exit 1
fi
if [ -z "$SAMPLE" ]; then # If no SAMPLE name specified then do not proceed
	usage
	echo "#ERROR: You need to specify a SAMPLE name that refers to your .bam file \$SAMPLE.marked.sort.bwa.$BUILD.bam."
	exit 1
fi
if [ -z "$WORKDIR" ]; then # If no output directory then use current directory
	WORKDIR=/hpcfs/users/${USER}/tmp/$SAMPLE
	echo "Using /hpcfs/users/${USER}/tmp/$SAMPLE as the output directory"
fi

tmpDir=/hpcfs/users/${USER}/tmp/$OUTPREFIX # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi
	
# Define files for the array
bedFile=($arrIndexBedFiles)

## Start of the script ##
# ApplyBQSR replaces PrintReads in GATK4 for the application of base quality score recalibration
# You should only run ApplyBQSR with the covariates table created from the input BAM
# -bqsr NOT -BQSR
# -RF OR --read-filter  NOT -rf BadCigar unrecognized read filter name
 
cd $tmpDir
java -Xmx6g -Djava.io.tmpdir=$tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]} -jar $GATKPATH/GenomeAnalysisTK.jar ApplyBQSR \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-I $BAMPATH/$SAMPLE.marked.sort.bwa.$BUILD.bam \
-L $ChrIndexPath/${bedFile[$SLURM_ARRAY_TASK_ID]} \
-bqsr $tmpDir/$SAMPLE.recal.grp \
--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
--emit-original-quals true \
-O ${bedFile[$SLURM_ARRAY_TASK_ID]}.$SAMPLE.realigned.recal.sorted.bwa.$BUILD.bam > $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$SAMPLE.pipeline.log 2>&1
