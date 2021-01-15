#!/bin/bash

#SBATCH -J MapSortDup
#SBATCH -o /hpcfs/users/%u/log/Dedup-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=06:00:00
#SBATCH --mem=124GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# run the executable
# A script to map reads and then call variants using the GATK v4.x best practices designed for the Phoenix supercomputer but will work on stand alone machines too

usage()
{
echo "# Script for processing and mapping Illumina 100bp pair-end sequence data and optionally plotting coverage for an interval
# Requires: BWA 0.7.x, Picard, samtools, GATKv4.x, BWA-Picard-GATK-CleanUp.sh.
# This script assumes your sequence files are gzipped
#
# Usage sbatch $0 -p file_prefix -s /path/to/sequences -o /path/to/output -c /path/to/config.cfg -S SAMPLE -L LIBRARY -I ID] | [ - h | --help ]
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
# Original: Derived from Illumina-Phred33-PE-FASTX-BWA-Picard-GATKv2.sh by Mark Corbett, 17/03/2014
# Modified: (Date; Name; Description)
# 24/09/2015; Mark Corbett; Fork original for genomes
# 25/09/2015; Mark Corbett; Pipe to samtools sort; Add getWGSMetrics
# 12/10/2015; Mark Corbett; Fix error collecting .bam files to merge
# 13/05/2016; Mark Corbett; Add option to specify sample name different from outPrefix.  Make seq file search explicit for *.fastq.gz
# 01/07/2016; Mark Corbett; Improve error handling
# 24/08/2016; Mark Corbett; Fork for HPC version, bring up to date with GATKv3.6
# 18/11/2016; Mark Corbett; Step down number of splits for PrintReads for higher efficiency
# 17/05/2017; Atma Ivancevic; translating for SLURM
# 16/11/2017; Mark Corbett; Fork to split scripts to see if this works on Phoenix; swap Picard for sambamba for duplicate marking, change workDir default
# 17/11/2020; Mark Corbett; Fork this script for mapping only
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
   config=/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/GATK4/BWA-GATKHC.GRCh37_full_analysis_set_phoenix.cfg
fi
source $config

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
if [ -z "$sample" ]; then # If sample name not specified then use "outPrefix"
	sample=$outPrefix
	echo "Using $outPrefix for sample name"
fi
if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=/hpcfs/users/${USER}/alignments/$sample
	echo "Using $workDir as the output directory"
fi
if [ -z "$LB" ]; then # If library not specified then use "IlluminaGenome"
	LB=IlluminaGenome
	echo "Using $LB for library name"
fi

tmpDir=/hpcfs/groups/phoenix-hpc-neurogenetics/tmp/BWA/$outPrefix # Use a tmp directory for all of the GATK and samtools temp files
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
seqFile1=$(ls *.fastq.gz | grep $outPrefix\_ | head -n 1) # Assume sequence files are some form of $outPrefix_fastq.gz
if [ -f "$seqFile1" ]; then
    seqFile2=$(ls *.fastq.gz | grep $outPrefix\_ | tail -n 1)
	fileCount=$(ls *.fastq.gz | grep $outPrefix\_ | wc -l | sed 's/[^0-9]*//g')
	if [ $fileCount -ne "2" ]; then
        echo "#WARN: I've found $fileCount sequence files but I was hoping for only 2. The R1 and R2 files will be concatenated before mapping, see below for details."
        cat $(ls *.fastq.gz | grep $outPrefix\_ | grep _R1) > $tmpDir/$outPrefix.cat_R1.fastq.gz
        seqPath=$tmpDir
		seqFile1=$outPrefix.cat_R1.fastq.gz
        echo "#INFO: The following R1 files were concatenated 
        $(ls *.fastq.gz | grep $outPrefix\_ | grep _R1)"
        cat $(ls *.fastq.gz | grep $outPrefix\_ | grep _R2) > $tmpDir/$outPrefix.cat_R2.fastq.gz
        seqFile2=$outPrefix.cat_R2.fastq.gz
        echo "#INFO: The following R2 files were concatenated 
        $(ls *.fastq.gz | grep $outPrefix\_ | grep _R2)"
	fi	
else
	seqFile1=$(ls *.fastq.gz | grep -w $outPrefix | head -n 1)
	seqFile2=$(ls *.fastq.gz | grep -w $outPrefix | tail -n 1)
	fileCount=$(ls *.fastq.gz | grep -w $outPrefix | wc -l | sed 's/[^0-9]*//g') # Otherwise try other seq file name options
	if [ $fileCount -ne "2" ]; then
        echo "#WARN: I've found $fileCount sequence files but I was hoping for only 2. The R1 and R2 files will be concatenated before mapping, see below for details."
        cat $(ls *.fastq.gz | grep  -w $outPrefix | grep _R1) > $tmpDir/$outPrefix.cat_R1.fastq.gz
        seqPath=$tmpDir
		seqFile1=$outPrefix.cat_R1.fastq.gz
        echo "#INFO: The following R1 files were concatenated 
        $(ls *.fastq.gz | grep  -w $outPrefix | grep _R1)"
        cat $(ls *.fastq.gz | grep  -w $outPrefix | grep _R2) > $tmpDir/$outPrefix.cat_R2.fastq.gz
        seqFile2=$outPrefix.cat_R2.fastq.gz
        echo "#INFO: The following R2 files were concatenated 
        $(ls *.fastq.gz | grep  -w $outPrefix | grep _R2)"
	fi
fi
if [ ! -f "$seqPath/$seqFile1" ]; then # Proceed to epic failure if can't locate unique seq file names
	echo "Sorry I can't find your sequence files! I'm using $outPrefix as part of the filename to locate them"
	exit 1
fi
if [ -z "$ID" ]; then
	ID=$(zcat $seqPath/$seqFile1 | head -n 1 | awk -F : '{OFS="."; print substr($1, 2, length($1)), $2, $3, $4}').$outPrefix # Hopefully unique identifier INSTRUMENT.RUN_ID.FLOWCELL.LANE.DNA_NUMBER. Information extracted from the fastq
fi

# load modules
module load arch/haswell
module load BWA/0.7.17-foss-2016b
module load HTSlib/1.10.2-foss-2016b
module load SAMtools/1.10-foss-2016b
module load sambamba/0.6.6-foss-2016b

## Start of the script ##
# Map reads to genome using BWA-MEM
# do NOT use the -M option when aligning to GRCh38 as alignment is alt-aware
# -K flag asks bwa-mem to load a fixed number of bases into RAM so enables reproducibility
 
cd $tmpDir
bwa mem -M -t 24 -R "@RG\tID:$ID\tLB:$LB\tPL:ILLUMINA\tSM:$sample" $BWAINDEXPATH/$BWAINDEX $seqPath/$seqFile1 $seqPath/$seqFile2 |\
samtools view -bT $GATKREFPATH/$GATKINDEX - |\
samtools sort -l 5 -m 4G -@24 -T$sample -o $sample.samsort.bwa.$BUILD.bam -

# Mark duplicates
sambamba markdup -t 24 -l 5 --tmpdir=$tmpDir --overflow-list-size 1000000 --hash-table-size 1000000 $sample.samsort.bwa.$BUILD.bam $workDir/$sample.marked.sort.bwa.$BUILD.bam
if [ -f "$workDir/$sample.marked.sort.bwa.$BUILD.bam" ]; then
    rm -R $tmpDir
else
	echo "Duplicate marking or earlier stage failed!"
	exit 1
fi

echo "# Flagstats" > $workDir/$sample.Stat_Summary.txt
samtools flagstat $workDir/$sample.marked.sort.bwa.$BUILD.bam >> $workDir/$sample.Stat_Summary.txt
