#!/bin/bash

#SBATCH -J STARmap
#SBATCH -o /hpcfs/users/%u/log/STARmap-slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1               	                                # number of nodes
#SBATCH -n 10              	                                # number of cores
#SBATCH --time=03:00:00    	                                # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=56G         	                                # memory pool for all cores

# Notification configuration 
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL   					# Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	# Email to which notification will be sent

# STAR.map.sh
#Set some paths and define functions
if [ -z "${enviroCfg}" ]; then # Test if the script was executed independently of the Universal Launcher script
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,rnaseq-germline-snps-indels,configs,g')"
    enviroCfg="${configDir}/BWA-GATKHC.environment.cfg"
fi
source ${enviroCfg}
threads=8 # Set one less than n above

usage()
{
echo "# STAR.map.sh.sh a slurm submission script for mapping Illumina paired end RNA-seq reads with STAR. 
# Before running as a batch script you need to know the number of samples you have.
#
# Dependencies:  An input text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2 /path/to/optional_SV_file\"
#                The SLURM log directory must exist ${userDir}/log or submission to SLURM will fail
#
# Usage: sbatch --array 0-(n Samples-1) $0 -i inputFile.txt [-o /path/to/outDir -c /path/to/config.cfg ] | $0 [-h | --help]
# If you have a lot of jobs you can limit the number running at a time to e.g. 8 like this:
# Usage: sbatch --array 0-(n Samples-1)%8 $0 -i inputFile.txt [ -o /path/to/outDir -c /path/to/config.cfg ] | $0 [-h | --help]
#
# Options: 
# -i	REQUIRED. Path and file name of a text file with sequences listed in the form \"read-group-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2 /path/to/optional_SV_file\"
#                 The optional SV file can be a VCF or a tab-delimited file as specified in the ARRIBA documents: https://arriba.readthedocs.io/en/latest/input-files/#structural-variant-calls-from-wgs
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
    # -i	REQUIRED. Path and file name of a text file with sequences listed in the form \"sample-ID path/to/read_1-1,...,path/to/read_n-1 /path/to/read_1-2,...,/path/to/read_n-2 /path/to/optional_SV_file\""
	exit 1
fi
if [ -z "${Config}" ]; then # If no Config file specified use the default
    Config=${scriptDir}/configs/GATK.RNAseq.germline.hg38.phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source ${Config}

# Define variables for the array jobs
sampleID=($(awk -F" " '{print $1}' ${SeqFile}))
read1=($(awk -F" " '{print $2}' ${SeqFile}))
read2=($(awk -F" " '{print $3}' ${SeqFile}))

if [ ! -d "${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}" ]; then
    mkdir -p ${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}
fi

tmpDir="${outDir}/${sampleID[$SLURM_ARRAY_TASK_ID]}/${sampleID[$SLURM_ARRAY_TASK_ID]}_STARtmp"
if [ -d "${tmpDir}" ]; then
    rm -r "${tmpDir}"  # STAR likes to start with a clean tmp directory
    echo "## WARN: A pre-existing STAR temporary directory ${tmpDir} was removed before starting this run.  This is usually occurs when a previous run failed."
fi

# Do the thing!
$STAR_prog \
    --runThreadN $threads \
    --genomeDir $STAR_index_dir/$buildID --genomeLoad NoSharedMemory --outTmpDir $tmpDir\
    --readFilesIn ${read1[$SLURM_ARRAY_TASK_ID]} ${read2[$SLURM_ARRAY_TASK_ID]} --readFilesCommand "gunzip -c" \
    --outSAMtype BAM SortedByCoordinate \
    ${"--sjdbOverhang "+(read_length-1)} \
    --twopassMode Basic \
    --limitBAMsortRAM ${star_mem+"000000000"} \
    --limitOutSJcollapsed ${default=1000000 star_limitOutSJcollapsed} \
    --outFileNamePrefix $outDir/${sampleID[$SLURM_ARRAY_TASK_ID]}/
