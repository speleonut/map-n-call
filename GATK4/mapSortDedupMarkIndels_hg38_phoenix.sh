#!/bin/bash

#SBATCH -J MapSortDup
#SBATCH -o /fast/users/%u/launch/slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=15:00:00
#SBATCH --mem=124GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# load modules
module load BWA/0.7.17-foss-2016b
module load Java/1.8.0_121
module load HTSlib/1.10.2-foss-2016b
module load SAMtools/1.10-foss-2016b
module load sambamba/0.6.6-foss-2016b
### module load GATK/4.0.0.0-Java-1.8.0_121
### module load picard/2.6.0 or higher

# run the executable
# A script to map reads and then call variants using the GATK v4.x best practices designed for the Phoenix supercomputer but will work on stand alone machines too

usage()
{
echo "# Script for processing and mapping Illumina 100bp pair-end sequence data and optionally plotting coverage for an interval
# Requires: BWA 0.7.x, Picard, samtools, sambamba, GATKv4.x, BWA-Picard-GATK-CleanUp.sh.
# This script assumes your sequence files are gzipped
#
# Usage $0 -p file_prefix -s /path/to/sequences -o /path/to/output -c /path/to/config.cfg -S SAMPLE -L LIBRARY -I ID] | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX_R1.fastq.gz 
# -s	REQUIRED. Path to the sequence files
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $FASTDIR/BWA-GATK/\$SAMPLE is used)
# -S	OPTIONAL. Sample name which will go into the BAM header and VCF header. If not specified, then will be set the same as -p
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line plain text, eg. MySeqProject20170317-PintOGuiness). Default \"IlluminaGenome\"
# -I	OPTIONAL. ID for the sequence (to go into the @RG line). If not specified the script will make one up from the first read header, and sample name
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
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
# 16/11/2017; Mark Corbett; Fork to split scripts to see if this works on Phoenix; swap Picard for sambamba for duplicate marking, change WORKDIR default
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
		-s )			shift
					SEQPATH=$1
					;;
		-S )			shift
					SAMPLE=$1
					;;
		-o )			shift
					WORKDIR=$1
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
if [ -z "$CONFIG" ]; then # If no config file specified use the default
   CONFIG=/data/neurogenetics/git/PhoenixScripts/shared/scripts/BWA-GATKHC.GRCh38_full_analysis_set_phoenix.cfg
fi
source $CONFIG

if [ -z "$OUTPREFIX" ]; then # If no file prefix specified then do not proceed
	usage
	echo "#ERROR: You need to specify a file prefix (PREFIX) referring to your sequence files eg. PREFIX_R1.fastq.gz."
	exit 1
fi
if [ -z "$SEQPATH" ]; then # If path to sequences not specified then do not proceed
	usage
	echo "#ERROR: You need to specify the path to your sequence files"
	exit 1
fi
if [ -z "$SAMPLE" ]; then # If sample name not specified then use "OUTPREFIX"
	SAMPLE=$OUTPREFIX
	echo "Using $OUTPREFIX for sample name"
fi
if [ -z "$WORKDIR" ]; then # If no output directory then use default directory
	WORKDIR=$FASTDIR/BWA-GATK/$SAMPLE
	echo "Using $FASTDIR/BWA-GATK/$SAMPLE as the output directory"
fi
if [ -z "$LB" ]; then # If library not specified then use "IlluminaGenome"
	LB=IlluminaGenome
	echo "Using $LB for library name"
fi

tmpDir=$FASTDIR/tmp/$OUTPREFIX # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi
	
# Locate sequence file names.
# This is a bit awkward and prone to errors since relies on only a few file naming conventions and assumes how they will line up after ls of files
# ...and assumes only your seq files are in the folder matching the file prefix
cd $SEQPATH
SEQFILE1=$(ls *.fastq.gz | grep $OUTPREFIX\_ | head -n 1) # Assume sequence files are some form of $OUTPREFIX_fastq.gz
if [ -f "$SEQFILE1" ]; then
	fileCount=$(ls *.fastq.gz | grep $OUTPREFIX\_ | wc -l | sed 's/[^0-9]*//g')
	if [ $fileCount -ne "2" ]; then
		echo "Sorry I've found the wrong number of sequence files and there's a risk I will map the wrong ones!"
		exit 1
	fi
	SEQFILE2=$(ls *.fastq.gz | grep $OUTPREFIX\_ | tail -n 1)
else
	fileCount=$(ls *.fastq.gz | grep -w $OUTPREFIX | wc -l | sed 's/[^0-9]*//g') # Otherwise try other seq file name options
	if [ $fileCount -ne "2" ]; then
		echo "Sorry I've found the wrong number of sequence files and there's a risk I will map the wrong ones!"
		exit 1
	fi
	SEQFILE1=$(ls *.fastq.gz | grep -w $OUTPREFIX | head -n 1) 
	SEQFILE2=$(ls *.fastq.gz | grep -w $OUTPREFIX | tail -n 1)
fi
if [ ! -f "$SEQFILE1" ]; then # Proceed to epic failure if can't locate unique seq file names
	echo "Sorry I can't find your sequence files! I'm using $OUTPREFIX as part of the filename to locate them"
	exit 1
fi
if [ -z $ID ]; then
	ID=$(zcat $SEQFILE1 | head -n 1 | awk -F : '{OFS="."; print substr($1, 2, length($1)), $2, $3, $4}').$OUTPREFIX # Hopefully unique identifier INSTRUMENT.RUN_ID.FLOWCELL.LANE.DNA_NUMBER. Information extracted from the fastq
fi

## Start of the script ##
# Map reads to genome using BWA-MEM
# do NOT use the -M option when aligning to GRCh38 as alignment is alt-aware
# -K flag asks bwa-mem to load a fixed number of bases into RAM so enables reproducibility
 
cd $tmpDir
bwa mem -K 100000000 -t 24 -R "@RG\tID:$ID\tLB:$LB\tPL:ILLUMINA\tSM:$SAMPLE" $BWAINDEXPATH/$BWAINDEX $SEQPATH/$SEQFILE1 $SEQPATH/$SEQFILE2 |\
samtools view -bT $GATKREFPATH/$GATKINDEX - |\
samtools sort -l 5 -m 4G -@24 -T$SAMPLE -o $SAMPLE.samsort.bwa.$BUILD.bam -

# Mark duplicates
sambamba markdup -t 24 -l 5 --tmpdir=$tmpDir --overflow-list-size 1000000 --hash-table-size 1000000 $SAMPLE.samsort.bwa.$BUILD.bam $SAMPLE.marked.sort.bwa.$BUILD.bam
if [ -f "$SAMPLE.marked.sort.bwa.$BUILD.bam" ]; then
    rm $SAMPLE.samsort.bwa.$BUILD.bam
else
	echo "Duplicate marking or earlier stage failed!"
	exit 1
fi

echo "# Flagstats" > $WORKDIR/$SAMPLE.Stat_Summary.txt
samtools flagstat $SAMPLE.marked.sort.bwa.$BUILD.bam >> $WORKDIR/$SAMPLE.Stat_Summary.txt

# Polish the BAM using GATK best practices
# IndelRealign NO longer included in GATK4, has been integrated into the variant callers
# Realign indels. eventually we should just output variants of interest with HC -bamOut option
# Create a file of intervals to realign to



# For the next steps split the bams into bits based on the IndexBedFiles
# First make tmp dirs
for bed in $arrIndexBedFiles; do
	mkdir -p $tmpDir/$bed
done
