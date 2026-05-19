#!/bin/bash
#SBATCH -J olDirtyB
#SBATCH -o /hpcfs/users/%u/log/GATK.RNA.germline.clean.up.slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1                            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 1                            # number of cores (here uses 12)
#SBATCH --time=00:10:00                 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=1G                        # memory pool for all cores (here set to 32 GB)

# Notification configuration
#SBATCH --mail-type=END                 # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  # Email to which notification will be sent

if [ -z "${enviroCfg}" ]; then # Test if the script was executed independently of the Universal Launcher script
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,rnaseq-germline-snps-indels,configs,g')"
    enviroCfg="${configDir}/BWA-GATKHC.environment.cfg"
fi
source ${enviroCfg}

usage()
{
echo "# clean.up.sh a slurm submission script for removing intermediate files from the GATK4 RNA germline variant calling pipeline. 
# Before running as a batch script you need to know the number of samples you have.
#
# Dependencies:  An input text file with one sampleID per line.  Other columns in the file will be ignored.
#                The SLURM log directory must exist ${userDir}/log or submission to SLURM will fail
#
# Usage: sbatch --array 0-(n Samples-1) $0 -i inputFile.txt [-o /path/to/outDir -c /path/to/config.cfg ] | $0 [-h | --help]
#
# Options: 
# -i	REQUIRED. Your original RNASeq sample input file. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -o	REQUIRED. Path to the root of the output directory which contains subfolders for each sampleID.
# -c	OPTIONAL. Path to a config file for the genome to be mapped. Default is GATK.RNAseq.germline.hg38.phoenix.cfg. 
# -h | --help     Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 19/05/2026
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
        -c )    shift
                Config=$1
                ;;
        -o )    shift
                outDir=$1
                ;;
        -h | --help )    usage
                         exit 0
                         ;;
        * )    usage
               exit 1
    esac
    shift
done

# Check that your script has everything it needs to start.
if [ -z "${seqFile}" ]; then #If sequence file list in a text file is not supplied then do not proceed
	usage
	echo "# ERROR: You need to specify the path and name of the sequence file list
    # -i	REQUIRED. Path and file name of a text file with one sample ID per line."
	exit 1
fi
if [ -z "${outDir}" ]; then #If the output directory is not supplied then do not proceed
	usage
	echo "# ERROR: You need to specify the path and name of the output directory
    # -o	REQUIRED. Path to the root of the output directory which contains subfolders for each sampleID."
	exit 1
fi

if [ -z "${Config}" ]; then # If no Config file specified use the default
    Config=${scriptDir}/configs/GATK.RNAseq.germline.hg38.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source ${Config}

# Define variables for the array jobs
sampleID=($(awk -F" " '{print $1}' ${seqFile}))

echo "## INFO: Removing the following intermediate files were removed.
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.${STARINDEX}.recal.split.marked.sort.bai
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.${STARINDEX}.recal.split.marked.sort.bam
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.marked.sort.bai
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.marked.sort.bam
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.split.marked.sort.bai
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.split.marked.sort.bam
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}_tmp
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/Aligned.sortedByCoord.out.bam
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/_STARgenome
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/_STARpass1"

rm -rf ${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.${STARINDEX}.recal.split.marked.sort.bai \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.${STARINDEX}.recal.split.marked.sort.bam \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.marked.sort.bai \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.marked.sort.bam \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.split.marked.sort.bai \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.split.marked.sort.bam \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}_tmp \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/Aligned.sortedByCoord.out.bam \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/_STARgenome \
${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/_STARpass1"
