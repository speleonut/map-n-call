#!/bin/bash
# This is the master script that coordinates job submission for the neurogenetics BWA-GATK4 haplotype caller pipeline.
## Set hard-coded paths and define functions ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call"
logDir="/hpcfs/users/${USER}/log"

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
fi

test_genome_build() {
case "${BUILD}" in
    "hs38DH" | "GRCh38_full_analysis_set" )    genomeType="has_alt_contigs"
                ;;
    "hs37d5" | "ucsc.hg19" )    genomeType="no_alt_contigs"
                ;;
    "CHM13v2" )    genomeType="no_alt_contigs"
                ;;
    * )         genomeType="has_alt_contigs"
                Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
                echo "## WARN: Genome build ${BUILD} not recognized, the default pipeline for GRCh38 / hg38 will be used."
                ;;
esac
}

usage()
{
echo "# This is the master script that coordinates job submission for primarily Illumina genome sequencing alignments but will work for exomes too.
# The script will select the right parameters to work with either GRCh37/hg19 or GRCh38/hg38 genome builds.  
# The scripts deliver an indel realigned BAM file, a WGS metrics report, a pipeline log and a gzipped gVCF file.
# Requires: BWA-MEM, samtools, sambamba, GATK, Picard, Java 
#
# Usage $0 -p file_prefix -s /path/to/sequences -o /path/to/output -c /path/to/config.cfg -S Sample -L LIBRARY -I ID] | [ - h | --help ]
#
# Options
# -p	REQUIRED. A prefix to your sequence files of the form PREFIX_R1.gz 
# -s	REQUIRED. Path to the sequence files
# -c	OPTIONAL. /path/to/config.cfg. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory /hpcfs/users/${USER}/BWA-GATK/\${Sample} is used)
# -S	OPTIONAL. Sample name which will go into the BAM header and VCF header. If not specified, then will be set the same as -p but it is recommended to set this if you can
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line plain text, eg. MySeqProject20170317-PintOGuiness). Default \"IlluminaGenome\"
# -I	OPTIONAL. ID for the sequence (to go into the @RG line). If not specified the script will make one up from the first read header, and Sample name
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Mark Corbett, 16/11/2017 
# Modified: (Date; Name; Description)
# 20/4/2020; Ali Gardner; for hg38 Phoenix
# 
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

## Pre-flight checks ##
if [ -z "$Config" ]; then # If no config file specified use the default
    Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config
if [ -z "${outPrefix}" ]; then # If no file prefix specified then do not proceed
	usage
	echo "## ERROR: You need to specify a file prefix (PREFIX) referring to your sequence files eg. PREFIX_R1.fastq.gz."
	exit 1
fi
if [ -z "$seqPath" ]; then # If path to sequences not specified then do not proceed
	usage
	echo "## ERROR: You need to specify the path to your sequence files"
	exit 1
fi

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

if [ -z "$LB" ]; then # If library not specified then use "IlluminaGenome"
	LB=IlluminaGenome
	echo "## INFO: Using $LB for library name"
fi
if [ -z "$Sample" ]; then # If Sample name not specified then use "outPrefix"
	Sample=${outPrefix}
	echo "## INFO: Using ${outPrefix} for Sample name"
fi
if [ -z "$workDir" ]; then # If no output directory then set and create a default directory
	workDir=/hpcfs/users/${USER}/BWA-GATKHC/$Sample
	echo "## INFO: Using $workDir as the output directory"
fi
if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi

## Launch the job chain ##
test_genome_build
case "${genomeType}" in
    "has_alt_contigs" ) BWAjob=`sbatch --export=ALL $scriptDir/GATK4/mapSortDedupMarkIndels_alt_aware_phoenix.sh  -c $Config -p ${outPrefix} -s $seqPath -S $Sample -o $workDir -L $LB -I $ID`
        ;;
    "no_alt_contigs" ) BWAjob=`sbatch --export=ALL $scriptDir/GATK4/mapSortDedupMarkIndels_no_alt_phoenix.sh  -c $Config -p ${outPrefix} -s $seqPath -S $Sample -o $workDir -L $LB -I $ID`
        ;;
    * ) echo "## ERROR: Genome build ${BUILD} corresponds with unrecognised pipeline version ${genomeType} the following config file is in use ${Config}"
        exit 1
esac
BWAjob=$(echo $BWAjob | cut -d" " -f4)
BQSRjob=`sbatch --export=ALL --dependency=afterok:${BWAjob} $scriptDir/GATK4/GATK.BQSR_Universal_phoenix.sh -c $Config -S $Sample -o $workDir`
BQSRjob=$(echo $BQSRjob | cut -d" " -f4)
ApplyBQSRJob=`sbatch --array=0-23 --export=ALL --dependency=afterok:${BQSRjob} $scriptDir/GATK4/GATK.ApplyBQSR_Universal_phoenix.sh -c $Config -S $Sample -o $workDir`
ApplyBQSRJob=$(echo $ApplyBQSRJob | cut -d" " -f4)
MergeJob=`sbatch --export=ALL --dependency=afterok:${ApplyBQSRJob} $scriptDir/GATK4/sambamba.Merge_Universal_phoenix.sh -c $Config -S $Sample -o $workDir`
MergeJob=$(echo $MergeJob | cut -d" " -f4)
GATKHCjob=`sbatch --array=0-23 --export=ALL --dependency=afterok:${MergeJob} $scriptDir/GATK4/GATK.HC_Universal_phoenix.sh -c $Config -S $Sample -o $workDir`
GATKHCjob=$(echo $GATKHCjob | cut -d" " -f4)
sbatch --export=ALL --dependency=afterok:${MergeJob} $scriptDir/GATK4/PicardCollectWGSMetrics_Universal_phoenix.sh -c $Config -S $Sample -o $workDir
sbatch --export=ALL --dependency=afterok:${GATKHCjob} $scriptDir/GATK4/GATK.gatherVCFs_Universal_phoenix.sh -c $Config -S $Sample -o $workDir
