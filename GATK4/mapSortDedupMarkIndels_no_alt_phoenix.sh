#!/bin/bash

#SBATCH -J MapSortDup
#SBATCH -o /hpcfs/users/%u/log/mapSortDedup-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 25
#SBATCH --time=15:00:00
#SBATCH --mem=148GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=${USER}@adelaide.edu.au

# run the executable
# A script to map reads and then call variants using the GATK v4.x best practices designed for the Phoenix supercomputer but will work on stand alone machines too
scriptDir=/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call
logDir="/hpcfs/users/${USER}/log"

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
fi

modList=("arch/haswell" "BWA/0.7.17-foss-2016b" "HTSlib/1.10.2-foss-2016b" "SAMtools/1.10-foss-2016b" "sambamba/0.6.6-foss-2016b")

usage()
{
echo "# Script for processing and mapping Illumina 100bp pair-end sequence data and optionally plotting coverage for an interval
# Requires: BWA 0.7.x, Picard, samtools, GATKv4.x, BWA-Picard-GATK-CleanUp.sh.
# This script assumes your sequence files are gzipped
#
# Usage: 
# sbatch $0 -p file_prefix -s /path/to/sequences -o /path/to/output -c /path/to/Config.cfg -S Sample -L LIBRARY -I ID] | [ - h | --help ]
# sbatch --array 0-(n-1Samples) $0 -i /path/to/multiSampleInputFile.txt -s /path/to/sequences -o /path/to/output -c /path/to/Config.cfg -L LIBRARY | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX_R1.fastq.gz 
# -s	REQUIRED. Path to the sequence files
# -i    OPTIONAL. An input file specifying prefix and Sample name in tab separated columns
# -c	OPTIONAL. /path/to/Config.cfg. Sets environment to match each genome build. A default Config will be used if this is not specified.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified, the output directory /hpcfs/users/${USER}/BWA-GATK/\$Sample is used)
# -S	OPTIONAL. Sample name which will go into the BAM header and VCF header. If not specified, then this will be the same as -p
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line plain text, eg. MySeqProject20170317-PintOGuiness). Default \"IlluminaLibrary\"
# -I	OPTIONAL. ID for the sequence (to go into the @RG line). If not specified the script will make one up from the first read header and the Sample name.
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
# 17/11/2020; Mark Corbett; Fork this script for mapping only
# 15/01/2021; Mark Corbett; Add in array job option
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
		-i )			shift
					inputFile=$1
					;;
		-h | --help )		usage
					exit 0
					;;
		* )			usage
					exit 1
	esac
	shift
done

# Check script requirements
if [ -z "$Config" ]; then # If no Config file specified use the default
   Config=$scriptDir/configs/BWA-GATKHC.GRCh37_full_analysis_set_phoenix.cfg
fi
source $Config
if [ ! -z "$inputFile" ]; then
    outPrefix=($(awk '{print $1}' $inputFile))
    Sample=($(awk '{print $2}' $inputFile))
fi
if [ -z "$outPrefix" ]; then # If no file prefix specified then do not proceed
	usage
	echo "## ERROR: You need to specify a file prefix referring to your sequence files eg. PREFIX_R1.fastq.gz."
	exit 1
fi
if [ -z "$seqPath" ]; then # If path to sequences is not specified then do not proceed
	usage
	echo "## ERROR: You need to specify the path to your sequence files"
	exit 1
fi
if [ -z "$Sample" ]; then # If Sample name not specified then use "outPrefix"
	Sample=${outPrefix[$SLURM_ARRAY_TASK_ID]}
	echo "## INFO: Using ${Sample} for Sample name"
fi
if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=/hpcfs/users/${USER}/alignments/${Sample[$SLURM_ARRAY_TASK_ID]}
	echo "#INFO: Using $workDir as the output directory"
fi
if [ -z "$LB" ]; then # If library not specified then use "IlluminaLibrary"
	LB=IlluminaLibrary
	echo "#INFO: Using $LB for library name"
fi

tmpDir=/hpcfs/groups/phoenix-hpc-neurogenetics/tmp/${USER}/${Sample[$SLURM_ARRAY_TASK_ID]} # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi
if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi	
# Locate sequence file names.
# This is a bit awkward and prone to errors since relies on only a few file naming conventions and assumes how they will line up after ls of files
# ...and assumes only your seq files are in the folder matching the file prefix
cd $seqPath
seqFile1=$(ls *.gz | grep ${outPrefix[$SLURM_ARRAY_TASK_ID]}\_ | head -n 1) # Assume sequence files are some form of $outPrefix_fastq.gz
if [ -f "$seqFile1" ]; then
    seqFile2=$(ls *.gz | grep ${outPrefix[$SLURM_ARRAY_TASK_ID]}\_ | tail -n 1)
	fileCount=$(ls *.gz | grep ${outPrefix[$SLURM_ARRAY_TASK_ID]}\_ | wc -l | sed 's/[^0-9]*//g')
	if [ $fileCount -ne "2" ]; then
        echo "## WARN: I've found $fileCount sequence files but I was hoping for only 2. The R1 and R2 files will be concatenated before mapping, see below for details."
        cat $(ls *.gz | grep ${outPrefix[$SLURM_ARRAY_TASK_ID]}\_ | grep _R1) > $tmpDir/${outPrefix[$SLURM_ARRAY_TASK_ID]}.cat_R1.fastq.gz
        seqPath=$tmpDir
		seqFile1=${outPrefix[$SLURM_ARRAY_TASK_ID]}.cat_R1.fastq.gz
        echo "## INFO: The following R1 files were concatenated 
        $(ls *.gz | grep ${outPrefix[$SLURM_ARRAY_TASK_ID]}\_ | grep _R1)"
        cat $(ls *.gz | grep ${outPrefix[$SLURM_ARRAY_TASK_ID]}\_ | grep _R2) > $tmpDir/${outPrefix[$SLURM_ARRAY_TASK_ID]}.cat_R2.fastq.gz
        seqFile2=${outPrefix[$SLURM_ARRAY_TASK_ID]}.cat_R2.fastq.gz
        echo "## INFO: The following R2 files were concatenated 
        $(ls *.gz | grep ${outPrefix[$SLURM_ARRAY_TASK_ID]}\_ | grep _R2)"
	fi	
else
	seqFile1=$(ls *.gz | grep -w ${outPrefix[$SLURM_ARRAY_TASK_ID]} | head -n 1)
	seqFile2=$(ls *.gz | grep -w ${outPrefix[$SLURM_ARRAY_TASK_ID]} | tail -n 1)
	fileCount=$(ls *.gz | grep -w ${outPrefix[$SLURM_ARRAY_TASK_ID]} | wc -l | sed 's/[^0-9]*//g') # Otherwise try other seq file name options
	if [ $fileCount -ne "2" ]; then
        echo "## WARN: I've found $fileCount sequence files but I was hoping for only 2. The R1 and R2 files will be concatenated before mapping, see below for details."
        cat $(ls *.gz | grep  -w ${outPrefix[$SLURM_ARRAY_TASK_ID]} | grep _R1) > $tmpDir/${outPrefix[$SLURM_ARRAY_TASK_ID]}.cat_R1.fastq.gz
        seqPath=$tmpDir
		seqFile1=${outPrefix[$SLURM_ARRAY_TASK_ID]}.cat_R1.fastq.gz
        echo "## INFO: The following R1 files were concatenated 
        $(ls *.gz | grep  -w ${outPrefix[$SLURM_ARRAY_TASK_ID]} | grep _R1)"
        cat $(ls *.gz | grep  -w ${outPrefix[$SLURM_ARRAY_TASK_ID]} | grep _R2) > $tmpDir/${outPrefix[$SLURM_ARRAY_TASK_ID]}.cat_R2.fastq.gz
        seqFile2=${outPrefix[$SLURM_ARRAY_TASK_ID]}.cat_R2.fastq.gz
        echo "## INFO: The following R2 files were concatenated 
        $(ls *.gz | grep  -w ${outPrefix[$SLURM_ARRAY_TASK_ID]} | grep _R2)"
	fi
fi
if [ ! -f "$seqPath/$seqFile1" ]; then # Proceed to epic failure if can't locate unique seq file names
	echo "## ERROR: Sorry I can't find your sequence files! I'm using ${outPrefix[$SLURM_ARRAY_TASK_ID]} as part of the filename to locate them"
	exit 1
fi
if [ -z "$ID" ]; then
	ID=$(zcat $seqPath/$seqFile1 | head -n 1 | awk -F : '{OFS="."; print substr($1, 2, length($1)), $2, $3, $4}').${outPrefix[$SLURM_ARRAY_TASK_ID]} # Hopefully unique identifier INSTRUMENT.RUN_ID.FLOWCELL.LANE.DNA_NUMBER. Information extracted from the fastq
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Start of the script ##
# Map reads to genome using BWA-MEM
 
cd $tmpDir
bwa mem  -K 100000000 -M -t 24 -R "@RG\tID:$ID\tLB:$LB\tPL:ILLUMINA\tSM:${Sample[$SLURM_ARRAY_TASK_ID]}" $BWAINDEXPATH/$BWAINDEX $seqPath/$seqFile1 $seqPath/$seqFile2 |\
samtools view -bT $GATKREFPATH/$BUILD/$GATKINDEX - |\
samtools sort -l 5 -m 4G -@24 -T${Sample[$SLURM_ARRAY_TASK_ID]} -o ${Sample[$SLURM_ARRAY_TASK_ID]}.samsort.bwa.$BUILD.bam -

# Mark duplicates 
sambamba markdup -t 24 -l 5 --tmpdir=$tmpDir --overflow-list-size 1000000 --hash-table-size 1000000 ${Sample[$SLURM_ARRAY_TASK_ID]}.samsort.bwa.$BUILD.bam $workDir/${Sample[$SLURM_ARRAY_TASK_ID]}.marked.sort.bwa.$BUILD.bam
if [ -f "$workDir/${Sample[$SLURM_ARRAY_TASK_ID]}.marked.sort.bwa.$BUILD.bam" ]; then
    rm -R $tmpDir
else
	echo "## ERROR: Duplicate marking or earlier stage failed!"
	exit 1
fi

echo "# Flagstats" > $workDir/${Sample[$SLURM_ARRAY_TASK_ID]}.Stat_Summary.txt
samtools flagstat $workDir/${Sample[$SLURM_ARRAY_TASK_ID]}.marked.sort.bwa.$BUILD.bam >> $workDir/${Sample[$SLURM_ARRAY_TASK_ID]}.Stat_Summary.txt
