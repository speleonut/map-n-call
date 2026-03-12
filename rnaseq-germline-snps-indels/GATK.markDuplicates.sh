#!/bin/bash

#SBATCH -J SupaDupa
#SBATCH -o /hpcfs/users/%u/log/SupaDupa-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1               	                                # number of nodes
#SBATCH -n 2              	                                # number of cores
#SBATCH --time=10:00:00    	                                # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=8G         	                                # memory pool for all cores

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	# Email to which notification will be sent

# GATK.markDuplicates.sh
#Set some paths and define functions
if [ -z "${enviroCfg}" ]; then # Test if the script was executed independently of the Universal Launcher script
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,rnaseq-germline-snps-indels,configs,g')"
    enviroCfg="${configDir}/BWA-GATKHC.environment.cfg"
fi
source ${enviroCfg}
threads=8 # Set one less than n above

modList=("Java/21.0.2" "Python/3.11.3-GCCcore-12.3.0")

usage()
{
echo "# GATK.markDuplicates.sh a slurm submission script for marking duplicates in RNA-seq BAM files. 
# Before running as a batch script you need to know the number of samples you have.
#
# Dependencies:  An input text file with one sampleID per line.  Other columns in the file will be ignored.
#                The SLURM log directory must exist ${userDir}/log or submission to SLURM will fail
#
# Usage: sbatch --array 0-(n Samples-1) $0 -i inputFile.txt [-o /path/to/outDir -c /path/to/config.cfg ] | $0 [-h | --help]
# If you have a lot of jobs you can limit the number running at a time to e.g. 8 like this:
# Usage: sbatch --array 0-(n Samples-1)%8 $0 -i inputFile.txt [ -o /path/to/outDir -c /path/to/config.cfg ] | $0 [-h | --help]
#
# Options: 
# -i	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2\"
# -c	OPTIONAL. Path to a config file for the genome to be mapped. Default is GATK.RNAseq.germline.hg38.phoenix.cfg. 
# -o	OPTIONAL. Path to where you want to find your files, the default is $userDir/RNASeq/arriba/genomeBuild/. 
#                 Each analyses will be put in a subfolder of this output directory using the sampleID if you use the default, or specifiy the directory yourself.
# -h | --help     Prints the message you are reading.
#
# History: 
# Script created by: Mark Corbett on 11/04/2024
# email: mark dot corbett is at adelaide.edu.au
# Modified (Date; Name; Description):
#
" 
}

# Parse script options
while [ "$1" != "" ]; do
    case $1 in
        -i )    shift
                SeqFile=$1
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

# Do the thing!
$GATKPATH/gatk --java-options "-Xmx8g -Djava.io.tmpdir=$tmpDir" \
    MarkDuplicates \
    --INPUT $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/Aligned.out.bam \
    --OUTPUT $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.marked.sort.bam  \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT \
    --METRICS_FILE $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.metrics \
    >> ${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}.${BUILD}.RNA.germline.pipeline.log 2>&1
