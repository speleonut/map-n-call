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


modList=("Java/21.0.2" "Python/3.11.3-GCCcore-12.3.0")

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
        -i )            shift
                        seqFile=$1
                        ;;
        -o )            shift
                        outDir=$1
                        ;;
        -h | --help )   usage
                        exit 0
                        ;;
        * )	            usage
                        exit 1
    esac
    shift
done

# Check that your script has everything it needs to start.
if [ -z "${SeqFile}" ]; then #If sequence file list in a text file is not supplied then do not proceed
	usage
	echo "# ERROR: You need to specify the path and name of the sequence file list
    # -i	REQUIRED. Path and file name of a text file with one sample ID per line."
	exit 1
fi
if [ -z "${Config}" ]; then # If no Config file specified use the default
    Config=${scriptDir}/configs/GATK.RNAseq.germline.hg38.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source ${Config}

# Define variables for the array jobs
sampleID=($(awk -F" " '{print $1}' ${SeqFile}))

if [ ! -d "${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}" ]; then
    mkdir -p ${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}
fi

tmpDir="${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}_tmp"
if [ ! -d "${tmpDir}" ]; then
    mkdir -p "${tmpDir}"
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

# Base quality score recalibration
$GATKPATH/gatk --java-options "-Xmx28g -Djava.io.tmpdir=$tmpDir" BaseRecalibrator \
-R ${GATKREFPATH}/${BUILD}/${GATKINDEX} \
-I ${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.split.marked.sort.bam \
--use-original-qualities \
--known-sites ${GATKREFPATH}/${BUILD}/${DBSNP} \
--known-sites ${GATKREFPATH}/${BUILD}/${OneKg_INDELS} \
--known-sites ${GATKREFPATH}/${BUILD}/${Mills_INDELS} \
--output $tmpDir/${sampleID[$SLURM_ARRAY_TASK_ID]}.recal.grp \
 >> ${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.${BUILD}.RNA.germline.pipeline.log 2>&1

echo "
# BQSR metrics" >> ${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.Stat_Summary.txt
cat ${tmpDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}.recal.grp >> ${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.Stat_Summary.txt
