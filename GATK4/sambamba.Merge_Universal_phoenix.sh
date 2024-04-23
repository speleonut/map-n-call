#!/bin/bash

#SBATCH -J Merge
#SBATCH -o /hpcfs/users/%u/log/mergeBAM-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=03:00:00
#SBATCH --mem=96GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# A script to merge bam files of the same Sample from multiple genomic intervals
## List modules and file paths ##
source ${enviroCfg}

usage()
{
echo "# A script to merge bam files of the same Sample from multiple genomic intervals
# Requires: sambamba
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
if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=${userDir}/alignments/$Sample
	echo "## INFO: Using $workDir as the output directory"
fi

# Merge with sambamba
find $tmpDir/*.$Sample.recal.sorted.bwa.$BUILD.bam > $tmpDir/$Sample.inputBAM.txt

$sambambaProg merge -t 24 -l 5 $workDir/$Sample.recal.sorted.bwa.$BUILD.bam $(cat $tmpDir/$Sample.inputBAM.txt)

# Clean up temporary .bam files
if [ -f "$workDir/$Sample.recal.sorted.bwa.$BUILD.bam" ]; then
	rm $tmpDir/*.$Sample.recal.sorted.bwa.$BUILD.ba* $workDir/$Sample.marked.sort.bwa.$BUILD.ba* 
else
    echo "## ERROR: Merging bam files failed please check log files for errors"
    exit 1
fi
