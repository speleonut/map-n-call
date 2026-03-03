#!/bin/bash

# A script to launch GATK4 Universal RNAseq Germline SNPs and Indels pipeline on Phoenix
## Hard coded paths for your system should be set in configs/BWA-GATKHC.environment.cfg ##
whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
configDir="$(echo ${whereAmI} | sed -e 's,rnaseq-germline-snps-indels,configs,g')"
utilitiesDir="$(echo ${whereAmI} | sed -e 's,rnaseq-germline-snps-indels,utilities,g')"
enviroCfg="${configDir}/BWA-GATKHC.environment.cfg"
source ${enviroCfg}

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
fi

usage() {
echo "# This script coordinates the launch of the GATK4 Universal RNAseq Germline SNPs and Indels pipeline on Phoenix. It requires an output directory and a config file as input. The config file should contain all necessary parameters for the pipeline, including paths to reference genomes, input data, and any specific settings for the GATK tools.
# Requirements: RNAseq fastq files or mapped (bam / cram) files.
#
# Usage: $0 -i /path/to/input.file.txt -o /path/to/output/directory/ -c /path/to/config/file.cfg | [-h | --help]
#
# Options:
# -i             REQUIRED: Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
#                          OR path to a text file with a list of paths to bam or cram files. Note this is slower to run.
#                          Don't mix the two formats in the same file. If you have both, split them into two files and run the pipeline separately for each file.
# -o             OPTIONAL: Output directory where the results will be stored. This directory will be created if it does not exist.
# -c             OPTIONAL: Path to the configuration file that contains all necessary parameters for the pipeline.
# -h | --help    OPTIONAL: Display this help message and exit.
#
# History: 
# Script created by: Mark Corbett on 02/03/2026
# email: mark dot corbett is at adelaide.edu.au
# Modified (Date; Name; Description):
#
"
}

# Parse script options
while [ "$1" != "" ]; do
	case $1 in
        -i )    shift
                seqFile=$1
			    ;;
        -o )    shift
                outDir=$1
                ;;
        -c )    shift
                config=$1
                ;;
        -h | --help )    usage
                         exit 0
                         ;;
        * )    usage
               exit 1
    esac
    shift
done

# Check the input file
if [ -z "${seqFile}" ]; then
    echo "## ERROR: No input file provided. Please provide a text file with the required format using the -i option."
    usage
    exit 1
fi  
numSamples=$(wc -l ${seqFile} | cut -d" " -f1)
numTasks=$((${numSamples}-1))

# Check the input file format and trust the user to have put the correct information in there (for now).
case $(awk -F"\t" 'NR==1{print NF}' ${seqFile}) in
    3)  echo "## INFO: Input file format suggests fastq file input."
        bamInput=false
        ;;
    1)  if [ -f $(head -n 1 ${seqFile}) ]; then
            bamInput=true
        else
            echo "## ERROR: The input file format suggests there should be a bam file but it does not exist.
            BAM_FILE_CHECKED:$(head -n 1 ${seqFile})"
            usage
            exit 1
        fi
        ;;
    *) echo "## ERROR: The input file format is not correct. Please provide a text file with the required format using the -i option."
       usage
       exit 1
       ;;
esac
if [ $(awk -F"\t" 'NR==1{print NF}' ${seqFile}) -ne 3 ]; then
    
fi

# Check the output directory and if it wasn't provided create a default directory.
if [ -z "${outDir}" ]; then
    outDir="${useDir}/variants/GATK4_universal_RNAseq_germline_snps_indels/$(date +%Y%m%d_%H%M%S)"
    echo "## INFO: No output directory provided. A default directory will be created at 
    OUTPUT_DIR:${outDir}."
fi

if [ ! -d "${outDir}" ]; then
    mkdir -p ${outDir}
    echo "## INFO: New output directory created, you'll find all of the results from this pipeline here:
    OUTPUT_DIR:${outDir}"
fi

# Check the config and if not specified use the default.
if [ -z "${config}" ]; then
    config="${configDir}/GATK4.RNAseq.germline.hg38.phoenix.cfg"
    echo "## INFO: No config file provided. A default config file will be used.
    CONFIG_FILE:${config}."
fi

# Start coordinating slum jobs
if [ "${bamInput}" = true ]; then
    touch ${outDir}/sequences/STAR.input.list.txt
    bam2fqJob=`sbatch --array=0-${numTasks} --export=ALL,enviroCfg=${enviroCfg} ${utilitiesDir}/bam2fq.samtools.sh -i ${seqFile} -o ${outDir}/sequences`
    bam2fqJobID=$(echo $bam2fqJob | cut -d" " -f4)
    makeSTARinputJob=`sbatch --dependency=afterok:${bam2fqJobID} --export=ALL,enviroCfg=${enviroCfg} ${whereAmI}/make.STAR.input.sh ${outDir}/sequences`
    makeSTARinputJobID=$(echo $makeSTARinputJob | cut -d" "-f4)
    seqFile=${outDir}/sequences/STAR.input.list.txt
    STARmapJob=`sbatch --array=0-${numTasks} --dependency=afterok:${makeSTARinputJobID} --export=ALL ${whereAmI}/STAR.map.sh -i ${seqFile} -o ${outDir} -c ${config}`
    STARmapJobID=$(echo $STARmapJob | cut -d" " -f4)
else
    STARmapJob=`sbatch --array=0-${numTasks} --export=ALL,enviroCfg=${enviroCfg} ${whereAmI}/STAR.map.sh -i ${seqFile} -o ${outDir} -c ${config}`
    STARmapJobID=$(echo $STARmapJob | cut -d" " -f4)
fi

# Steps remaining from 
# https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl
