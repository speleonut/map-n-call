#!/bin/bash

#SBATCH -J WGSMetrics
#SBATCH -o /fast/users/%u/launch/slurm-%j.out

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
module load Java/1.8.0_121
### module load GATK 4.1.6.0
### module load picard/2.6.0 or higher

# run the executable
# A script to map reads and then call variants using the GATK v4.x best practices designed for the Phoenix supercomputer but will work on stand alone machines too

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
   CONFIG=/data/neurogenetics/git/PhoenixScripts/shared/scripts/BWA-GATKHC.GRCh38_full_analysis_set_phoenix.cfg
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
	WORKDIR=$FASTDIR/BWA-GATK/$SAMPLE
	echo "Using $FASTDIR/BWA-GATK/$SAMPLE as the output directory"
fi

tmpDir=$FASTDIR/tmp/$OUTPREFIX # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi
 
cd $WORKDIR
java -Xmx8g -Djava.io.tmpdir=$tmpDir -jar $PICARDPATH/picard.jar CollectWgsMetrics \
INPUT=$WORKDIR/$SAMPLE.realigned.recal.sorted.bwa.$BUILD.bam \
OUTPUT=$WORKDIR/$SAMPLE.WGS.Metrics \
REFERENCE_SEQUENCE=$GATKREFPATH/$GATKINDEX \
COVERAGE_CAP=100 \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=12000000 > $WORKDIR/$SAMPLE.WGS.Metrics.pipeline.log 2>&1
