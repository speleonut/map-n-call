#!/bin/bash

#SBATCH -J MapSortDup
#SBATCH -o /hpcfs/users/%u/log/mapSortDedup-slurm-%j.out
#SBATCH -p skylake,icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 25
#SBATCH --time=24:00:00
#SBATCH --mem=148GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Script for mapping Illumina pair-end sequence data
## Set hard coded paths and variables ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call"
logDir="/hpcfs/users/${USER}/log"
sambambaProg=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/sambamba-0.8.2-linux-amd64-static

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
fi
module purge
module use /apps/skl/modules/all
modList=("BWA/0.7.17-GCCcore-11.2.0" "HTSlib/1.17-GCC-11.2.0" "SAMtools/1.17-GCC-11.2.0")
usage()
{
echo "# Script for mapping Illumina pair-end sequence data
# Requires: BWA 0.7.x, samtools, sambamba
# This script assumes your sequence files are gzipped
#
# Usage $0 -p file_prefix -s /path/to/sequences -o /path/to/output -c /path/to/config.cfg -S Sample -L LIBRARY -I ID] | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX_R1.gz 
# -s	REQUIRED. Path to the sequence files
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory /hpcfs/users/${USER}/BWA-GATK/\$Sample is used)
# -S	OPTIONAL. Sample name which will go into the BAM header and VCF header. If not specified, then will be set the same as -p
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line plain text, eg. MySeqProject20170317-PintOGuiness). Default \"IlluminaGenome\"
# -I	OPTIONAL. ID for the sequence (to go into the @RG line). If not specified the script will make one up from the first read header, and Sample name
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Derived from Illumina-Phred33-PE-FASTX-BWA-Picard-GATKv2.sh by Mark Corbett, 17/03/2014
# Modified: (Date; Name; Description)
# 24/09/2015; Mark Corbett; Fork original for genomes
# 25/09/2015; Mark Corbett; Pipe to samtools sort; Add getWGSMetrics
# 12/10/2015; Mark Corbett; Fix error collecting .bam files to merge
# 13/05/2016; Mark Corbett; Add option to specify Sample name different from outPrefix.  Make seq file search explicit for *.gz
# 01/07/2016; Mark Corbett; Improve error handling
# 24/08/2016; Mark Corbett; Fork for HPC version, bring up to date with GATKv3.6
# 18/11/2016; Mark Corbett; Step down number of splits for PrintReads for higher efficiency
# 17/05/2017; Atma Ivancevic; translating for SLURM
# 16/11/2017; Mark Corbett; Fork to split scripts to see if this works on Phoenix; swap Picard for sambamba for duplicate marking, change workDir default
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-c )			shift
					Config=$1
					;;
		-p )			shift
					outPrefix=$1
					;;
		-s )			shift
					seqPath=$1
					;;
		-S )			shift
					Sample=$1
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
					exit 0
					;;
		* )			usage
					exit 1
	esac
	shift
done
if [ -z "$Config" ]; then # If no Config file specified use the default
    Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config

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
if [ -z "$Sample" ]; then # If Sample name not specified then use "outPrefix"
	Sample=$outPrefix
	echo "##INFO: Using $outPrefix for Sample name"
fi
if [ -z "$workDir" ]; then # If no output directory then use default directory
    workDir=/hpcfs/users/${USER}/BWA-GATK/$Sample
fi

if [ -z "$LB" ]; then # If library not specified then use "IlluminaGenome"
	LB=IlluminaGenome
	echo "##INFO: Using $LB for library name"
fi

## Locate sequence files ##
# This is a bit awkward and prone to errors since relies on only a few file naming conventions and assumes that they will line up correctly after finding the files
# ...and assumes only your seq files are in the folder matching the file prefix
seqFile1=$(find ${seqPath}/*.gz | grep ${outPrefix}\_ | head -n 1) # Assume sequence files are some form of ${outPrefix}_*.gz
if [ -f "$seqFile1" ]; then
	fileCount=$(find ${seqPath}/*.gz | grep ${outPrefix}\_ | wc -l | sed 's/[^0-9]*//g')
	if [ $fileCount -ne "2" ]; then
		echo "Sorry I've found the wrong number of sequence files (${fileCount}) and there's a risk I will map the wrong ones!"
		exit 1
	fi
	seqFile2=$(find ${seqPath}/*.gz | grep ${outPrefix}\_ | tail -n 1)
else
	fileCount=$(find ${seqPath}/*.gz | grep -w ${outPrefix} | wc -l | sed 's/[^0-9]*//g') # Otherwise try other seq file name options
	if [ $fileCount -ne "2" ]; then
		echo "Sorry I've found the wrong number of sequence files (${fileCount}) and there's a risk I will map the wrong ones!"
		exit 1
	fi
	seqFile1=$(find ${seqPath}/*.gz | grep -w ${outPrefix} | head -n 1) 
	seqFile2=$(find ${seqPath}/*.gz | grep -w ${outPrefix} | tail -n 1)
fi
if [ ! -f "${seqFile1}" ]; then # Proceed to epic failure if can't locate unique seq file names
	echo "Sorry I can't find your sequence files! I'm using ${outPrefix} as part of the filename to locate them"
	exit 1
fi
if [ -z $ID ]; then
	ID=$(zcat $seqFile1 | head -n 1 | awk -F : '{OFS="."; print substr($1, 2, length($1)), $2, $3, $4}').$outPrefix # Hopefully unique identifier INSTRUMENT.RUN_ID.FLOWCELL.LANE.DNA_NUMBER. Information extracted from the fastq
fi

## Create essential directories ##
tmpDir=/hpcfs/groups/phoenix-hpc-neurogenetics/tmp/${USER}/${Sample} # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi

if [ ! -d "${workDir}"} ]; then
    mkdir -p $workDir
    echo "##INFO: Using $workDir as the output directory"
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Start of the script ##
# Map reads to genome using BWA-MEM
# do NOT use the -M option when aligning to GRCh38 as alignment is alt-aware
# -K flag asks bwa-mem to load a fixed number of bases into RAM so enables reproducibility
 
cd $tmpDir
bwa mem -K 100000000 -t 24 -R "@RG\tID:$ID\tLB:$LB\tPL:ILLUMINA\tSM:$Sample" $BWAINDEXPATH/$BWAINDEX $seqFile1 $seqFile2 |\
samtools view -bT $GATKREFPATH/$BUILD/$GATKINDEX - |\
samtools sort -l 5 -m 4G -@24 -T$Sample -o $Sample.samsort.bwa.$BUILD.bam -

# Mark duplicates
$sambambaProg markdup -t 24 -l 5 --tmpdir=$tmpDir --overflow-list-size 1000000 --hash-table-size 1000000 $Sample.samsort.bwa.$BUILD.bam $workDir/$Sample.marked.sort.bwa.$BUILD.bam
if [ -f "$workDir/$Sample.marked.sort.bwa.$BUILD.bam" ]; then
    rm $Sample.samsort.bwa.$BUILD.bam
else
	echo "## ERROR: Duplicate marking or earlier stage failed!"
	exit 1
fi

echo "# Flagstats" > $workDir/$Sample.Stat_Summary.txt
samtools flagstat $workDir/$Sample.marked.sort.bwa.$BUILD.bam >> $workDir/$Sample.Stat_Summary.txt
