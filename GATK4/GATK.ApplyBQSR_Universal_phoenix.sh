#!/bin/bash

#SBATCH -J ApplyBQSR
#SBATCH -o /hpcfs/users/%u/log/applyBQSR-slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=02:00:00
#SBATCH --mem=8GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# A script to apply base quality score recalibration using the GATK v4.x best practices

## List modules and file paths ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call"
modList=("arch/haswell" "Java/1.8.0_121")

usage()
{
echo "# A script to apply base quality score recalibration using the GATK v4.x best practices
# Requires: GATKv4.x
#
# Usage: sbatch --array 0-23 $0 -S sample_name [ -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
#
# Options
# -S	REQUIRED. Sample name if not specified then will be set the same as -p
# -c	OPTIONAL. /path/to/Config.cfg. A default Config will be used if this is not specified.  The Config contains all of the stuff that used to be set in the top part of our scripts
# -o	OPTIONAL. Path to where you want to find your file output (if not specified current directory is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# 
# Original: Derived from Illumina-Phred33-PE-FASTX-BWA-Picard-GATKv2.sh by Mark Corbett, 17/03/2014
# Modified: (Date; Name; Description)
# 24/09/2015; Mark Corbett; Fork original for genomes
# 25/09/2015; Mark Corbett; Pipe to samtools sort; Add getWGSMetrics
# 12/10/2015; Mark Corbett; Fix error collecting .bam files to merge
# 13/05/2016; Mark Corbett; Add option to specify Sample name different from OUTPREFIX.  Make seq file search explicit for *.fastq.gz
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
        -c )            shift
                        Config=$1
                        ;;
        -S )            shift
                        Sample=$1
                        ;;
        -o )            shift
                        workDir=$1
                        ;;
        -h | --help )   usage
                        exit 0
                        ;;
        * )	            usage
                        exit 1
    esac
    shift
done
if [ -z "$Config" ]; then # If no Config file specified use the default
    Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config
if [ -z "$Sample" ]; then # If no Sample name specified then do not proceed
	usage
	echo "## ERROR: You need to specify a Sample name that refers to your .bam file \$Sample.marked.sort.bwa.$BUILD.bam."
	exit 1
fi
if [ -z "$workDir" ]; then # If no output directory then use current directory
	workDir=/hpcfs/users/${USER}/BWA-GATK/$Sample
	echo "## INFO: Using $workDir as the output directory"
fi

tmpDir=/hpcfs/groups/phoenix-hpc-neurogenetics/tmp/${USER}/${Sample} # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi
	
## Define the array for the batch job ##
bedFile=($arrIndexBedFiles)

## Prepare for the next steps in the pipeline ##
# First make tmp dirs to hold files of each genome interval
for bed in $arrIndexBedFiles; do
	mkdir -p $tmpDir/$bed
done

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Start of the script ##
# ApplyBQSR replaces PrintReads in GATK4 for the application of base quality score recalibration
# You should only run ApplyBQSR with the covariates table created from the input BAM
 
cd $tmpDir
java -Xmx6g -Djava.io.tmpdir=$tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]} -jar $GATKPATH/GenomeAnalysisTK.jar ApplyBQSR \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-I $workDir/$Sample.marked.sort.bwa.$BUILD.bam \
-L $ChrIndexPath/${bedFile[$SLURM_ARRAY_TASK_ID]} \
-bqsr $tmpDir/$Sample.recal.grp \
--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
--emit-original-quals true \
-O ${bedFile[$SLURM_ARRAY_TASK_ID]}.$Sample.recal.sorted.bwa.$BUILD.bam >> $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$Sample.pipeline.log 2>&1
