#!/bin/bash

#SBATCH -J WGSMetrics
#SBATCH -o /hpcfs/users/%u/log/collectWGSmetrics-slurm-%j.out

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

# A script to collect alignment metrics from WGS bam files using Picard
## List modules and file paths ##
scriptDir="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/mark/map-n-call"
modList=("arch/haswell" "Java/1.8.0_121")

usage()
{
echo "# A script to collect alignment metrics from WGS bam files using Picard
# Requires: Picard, GATKv4.x
# This script assumes your sequence files are gzipped
#
# Usage $0 -S Sample_name [ -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
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

tmpDir=/hpcfs/users/${USER}/tmp/${Sample} # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Start of the script ##
cd $workDir
java -Xmx8g -Djava.io.tmpdir=$tmpDir -jar $PICARDPATH/picard.jar CollectWgsMetrics \
INPUT=$workDir/$Sample.recal.sorted.bwa.$BUILD.bam \
OUTPUT=$workDir/$Sample.$BUILD.WGS.Metrics \
REFERENCE_SEQUENCE=$GATKREFPATH/$BUILD/$GATKINDEX \
COVERAGE_CAP=100 \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM=12000000 > $workDir/$Sample.WGS.Metrics.pipeline.log 2>&1

rm -r $tmpDir
