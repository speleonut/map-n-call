#!/bin/bash

#SBATCH -J BQSR
#SBATCH -o /hpcfs/users/%u/log/bqsr-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=12:00:00
#SBATCH --mem=28GB

# Notification Configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# A script to calculate base quality score recalibrations using the GATK v4.x best practices

## List modules and file paths ##
source ${enviroCfg}


modList=("Java/17.0.6" "Python/3.11.3-GCCcore-12.3.0")

usage()
{
echo "# A script to calculate base quality score recalibrations using the GATK v4.x best practices
# Requires: GATKv4.x
#
# Usage $0 -S sample_name [ -o /path/to/output -c /path/to/config.cfg ] | [ - h | --help ]
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
# 13/05/2016; Mark Corbett; Add option to specify Sample name different from outPrefix.  Make seq file search explicit for *.fastq.gz
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

if [ -z "$Sample" ]; then # If no Sample name specified then do not proceed
	usage
	echo "## ERROR: You need to specify a Sample name that refers to your .bam file \$Sample.marked.sort.bwa.$BUILD.bam."
	exit 1
fi

if [ -z "${enviroCfg}" ]; then # Test if the script was executed independently of the Universal Launcher script
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,GATK4,configs,g')"
    source ${configDir}/BWA-GATKHC.environment.cfg
    if [ ! -d "${logDir}" ]; then
        mkdir -p ${logDir}
        echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
    fi
    tmpDir=${tmpDir}/${Sample}
    if [ ! -d "$tmpDir" ]; then
        mkdir -p $tmpDir
    fi
fi

if [ -z "$Config" ]; then # If no Config file specified use the default
    Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config

if [ -z "$workDir" ]; then # If no output directory then use current directory
	workDir=${userDir}/alignments/${Sample}
	echo "## INFO: Using $workDir as the output directory"
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Start data processing ##
cd $tmpDir

# Base quality score recalibration
$GATKPATH/gatk --java-options "-Xmx96g -Djava.io.tmpdir=$tmpDir" BaseRecalibrator \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-I $workDir/$Sample.marked.sort.bwa.$BUILD.bam \
--known-sites $GATKREFPATH/$BUILD/$DBSNP \
--known-sites $GATKREFPATH/$BUILD/${OneKg_INDELS} \
--known-sites $GATKREFPATH/$BUILD/${Mills_INDELS} \
--output $tmpDir/$Sample.recal.grp \
 >> $workDir/${Sample}.${BUILD}.pipeline.log 2>&1

echo "
# BQSR metrics" >> $workDir/$Sample.Stat_Summary.txt
cat $tmpDir/$Sample.recal.grp >> $workDir/$Sample.Stat_Summary.txt
