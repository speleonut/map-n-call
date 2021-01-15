#!/bin/bash
#SBATCH -J BAM2fq
#SBATCH -o /hpcfs/users/%u/log/bam2fq.samtools.slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8                            # number of cores (here uses 8)
#SBATCH --time=05:30:00                 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=32G                       # memory pool for all cores (here set to 32 GB)

# Notification configuration
#SBATCH --mail-type=END                 # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  # Email to which notification will be sent

#Script Paths
userDir=/hpcfs/users/${USER}

# bam2fq.samtools.sh
usage()
{
echo "# bam2fq.samtools.sh Sort a BAM by read name then convert to gzipped fastq files.
# Dependencies:  samtools v1.9+
# Info: http://www.htslib.org/doc/samtools.html
#
# Usage: sbatch $0 -b /path/to/bam/folder -o /path/to/output/folder -S sampleID | [-h | --help]
#        sbatch --array <0-(n-1) samples> $0 -i /path/to/input-file -o /path/to/output/folder
#
# Options:
# -b <arg>           REQUIRED: Path to where your bam file is located
# -S <arg>           REQUIRED: ID of the sample which will form the first part of your fastq file names
# -o <arg>           OPTIONAL: Path to the output default: $userDir/sequences/sampleID/SLURM_JOB_ID
# -h | --help        Prints the message you are reading.
#
# Array job options:
# -i <arg>           REQUIRED: If running an array job a tab delimited text file with two columns listing the required arguments -b and -S from above (in that order).
# -o <arg>           OPTIONAL: Path to the output default: $userDir/sequences/sampleID/SLURM_JOB_ID
# -h | --help        Prints the message you are reading.
#
# History:
# Script created by: Mark Corbett on 15/04/2020
# email:mark dot corbett is at adelaide university
# Modified (Date; Name; Description):
# 15/01/2021; Mark; Update for new Phoenix and create array job option
#
"
}

# Parse script options
while [ "$1" != "" ]; do
    case $1 in
        -b )    shift
                bamDir=$1
                ;;
        -S )    shift
                sampleID=$1
                ;;
        -o )    shift
                outDir=$1
                ;;
        -i )    shift
                inputFile=$1
                ;;
        -h | --help )   module load arch/haswell
                        module load SAMtools
                        samtools sort
                        samtools fastq
                        module unload SAMtools
                        module unload arch/haswell
                        usage
                        exit 0
                        ;;
        * )     module load arch/haswell
                module load SAMtools
                samtools sort
                samtools fastq
                module unload SAMtools
                module unload arch/haswell
                usage
                exit 1
    esac
    shift
done

# Check that your script has everything it needs to start.
if [ ! -z "$inputFile" ]; then
    bamDir=($(awk '{print $1}' $inputFile))
    sampleID=($(awk '{print $2}' $inputFile))
fi	
if [ -z "$bamDir" ]; then # If bamFile not specified then do not proceed
    usage
    echo "#ERROR: You need to specify -b /path/to/bamfile
    # -b <arg>    REQUIRED: Path to where your bam file is located"
    exit 1
fi
if [ -z "$sampleID" ]; then # If sample not specified then do not proceed
    usage
    echo "#ERROR: You need to specify -S sampleID because I need this to make your file names
    # -S <arg>    ID of the sample which will form the first part of your fastq file names"
    exit 1
fi

bamFile=$( find ${bamDir[SLURM_ARRAY_TASK_ID]}/*.bam | grep ${sampleID[SLURM_ARRAY_TASK_ID]} )

if [ -z "$outDir" ]; then # If output directory not specified then make one up
    outDir=$userDir/sequences/${sampleID[SLURM_ARRAY_TASK_ID]}/$SLURM_JOB_ID
    echo "#INFO: You didn't specify an output directory so I'm going to put your files here.
    $outDir"
fi
if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi
tmpDir=$outDir/tmp.$SLURM_JOB_ID
if [ ! -d "$tmpDir" ]; then
    mkdir -p $tmpDir
fi

# Load modules
module load arch/haswell
module load SAMtools

# Revert BAMs to fastq
samtools sort -l 0 -m 4G -n -@8 -T$tmpDir ${bamFile} |\
samtools fastq -1 $outDir/${sampleID[SLURM_ARRAY_TASK_ID]}.reads_R1.fastq.gz -2 $outDir/${sampleID[SLURM_ARRAY_TASK_ID]}.reads_R2.fastq.gz -0 /dev/null -s $outDir/${sampleID[SLURM_ARRAY_TASK_ID]}.reads_U1.fastq.gz -n -@8 -

# Clean up and generate run stats
rm -r $tmpDir
rcstat $SLURM_JOB_ID
