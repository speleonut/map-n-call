#!/bin/bash
#SBATCH -J getHugeBro
#SBATCH -o /hpcfs/users/%u/log/gunzip.slurm-%j.out

#SBATCH -p skylake,icelake,a100cpu
#SBATCH -N 1                            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 1                            # number of cores (here uses 8)
#SBATCH --time=08:00:00                 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=4G                       # memory pool for all cores (here set to 4 GB)

# Notification configuration
#SBATCH --mail-type=END                 # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  # Email to which notification will be sent

usage()
{
echo "# bulkUnzip.sh useful for when other crappy software does not handle this for you.
#
# Usage: 
# single sample: sbatch $0 -i /path/to/input-file | [ -h | --help ]
# multiple samples: sbatch --array <0-(n-1) samples> $0 -i /path/to/input-file | [ -h | --help ]
#
# Options:
# -i <arg>           REQUIRED: List of files with full or relative (used cautiously) file paths.
# -h | --help        Prints the message you are reading.
#
# History:
# Script created by: Mark Corbett on 01/07/2022
# email:mark dot corbett is at adelaide university
# Modified (Date; Name; Description):
# 10/10/2023; Mark; Add to map-n-call and update packages for new HPC
#
"
}

# Parse script options
while [ "$1" != "" ]; do
    case $1 in
        -i )    shift
                inputFile=$1
                ;;
        -h | --help )   usage
                        exit 0
                        ;;
        * )     usage
                exit 1
    esac
    shift
done

# Check that your script has everything it needs to start.
if [ -z "$inputFile" ]; then
    usage
    echo "#ERROR: You need to specify -i /path/to/input-file"
    exit 1
fi

module load gzip/1.10-GCCcore-11.2.0
readarray -t fileToDecompress < $inputFile

# Skywalker: "I can't. It's too big."" Yoda: "Size matters not. Look at me. Judge me by my size, do you? Hmm? Hmm."
gunzip ${fileToDecompress[$SLURM_ARRAY_TASK_ID]}
