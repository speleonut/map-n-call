#!/bin/bash
#SBATCH -J BAM2CRAM
#SBATCH -o /hpcfs/users/%u/log/bam2cram.samtools.slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch                        # partition (this is the queue your job will be added to)
#SBATCH -N 1                            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 8                            # number of cores (here uses 8)
#SBATCH --time=05:30:00                 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=36G                       # memory pool for all cores (here set to 32 GB)

# Notification configuration
#SBATCH --mail-type=END                 # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  # Email to which notification will be sent

#Script Paths
userDir=/hpcfs/users/${USER}

modList=("arch/haswell" "SAMtools/1.10-foss-2016b")

# bam2cram.samtools.sh
usage()
{
echo "# bam2cram.samtools.sh convert a BAM to CRAM file.
# Dependencies:  samtools v1.9+
# Info: http://www.htslib.org/doc/samtools.html
#
# Usage: sbatch $0 -b /path/to/bam/folder -S sampleID -g /path/to/reference -o /path/to/output/folder | [-h | --help]
#        sbatch --array <0-(n-1) samples> $0 -i /path/to/input-file -o /path/to/output/folder
#
# Options:
# -b <arg>           REQUIRED: Path to where your bam file is located.
# -S <arg>           REQUIRED: ID of the sample which will be used to identify your bam file.
# -g <arg>           REQUIRED: Path to the original reference that your BAM file was mapped to (Hint: this info is usually in the BAM file header)
# -o <arg>           OPTIONAL: Path to the output default: $userDir/alignments/sampleID/SLURM_JOB_ID.
# -h | --help        Prints the message you are reading.
#
# Array job options:
# -i <arg>           REQUIRED: If running an array job a tab delimited text file with two columns listing the required arguments -b and -S from above (in that order).
# -o <arg>           OPTIONAL: Path to the output default: $userDir/alignments/sampleID/SLURM_JOB_ID
# -h | --help        Prints the message you are reading.
#
# History:
# Script created by: Mark Corbett on 12/10/2021
# email:mark dot corbett is at adelaide university
# Modified (Date; Name; Description):
# 
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
        -g )    shift
                genome=$1
                ;;
        -o )    shift
                outDir=$1
                ;;
        -i )    shift
                inputFile=$1
                ;;
        -h | --help )   for mod in "${modList[@]}"; do
                            module load $mod
                        done
                        samtools view
                        module unload ${modList[1]}
                        module unload ${modList[0]}
                        usage
                        exit 0
                        ;;
        * )     for mod in "${modList[@]}"; do
                    module load $mod
                done
                samtools view
                module unload ${modList[1]}
                module unload ${modList[0]}
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
if [ -z "$genome" ]; then # If genome not specified then do not proceed
    usage
    echo "#ERROR: You need to specify -g /path/to/genome because I can't make a CRAM without this
    # -g <arg>    REQUIRED: Path to the original reference that your BAM file was mapped to (Hint: this info is usually in the BAM file header)"
    exit 1
fi

bamFile=$( find ${bamDir[SLURM_ARRAY_TASK_ID]}/*.bam | grep ${sampleID[SLURM_ARRAY_TASK_ID]} )
baseBamFile=$( basename ${bamFile} .bam )

if [ -z "$outDir" ]; then # If output directory not specified then make one up
    outDir=$userDir/alignments/${sampleID[SLURM_ARRAY_TASK_ID]}/$SLURM_JOB_ID
    echo "#INFO: You didn't specify an output directory so I'm going to put your files here.
    $outDir"
fi
if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi

# Load modules
for mod in "${modList[@]}"; do
    module load $mod
done

# Convert BAMs to CRAMs
samtools view -T ${genome} -C -@8 ${bamFile} -O CRAM -o ${outDir}/${baseBamFile}.cram
samtools index ${outDir}/${baseBamFile}.cram

# Send a polite reminder to clean up the old BAM file
echo "## INFO: Thanks for converting your BAM to a CRAM.  Help save some space by removing the original BAM file (you can always convert it back or remap it if you ever need the BAM again)"
