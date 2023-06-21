#!/bin/bash
#SBATCH -J MAS
#SBATCH -o /hpcfs/users/%u/log/revertSam.slurm-%j.out
#SBATCH -p skylake,icelake,skylakehm,v100cpu
#SBATCH -N 1                                                # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 4                                                # number of cores (here uses 2)
#SBATCH --time=08:00:00                                     # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=16G                                            # memory pool for all cores (here set to 8 GB)

# Notification configuration
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                		    # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	    # Email to which notification will be sent

# revertSam.sh
# Set location of picard.jar
PICARD=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/gatk-latest/GenomeAnalysisTK.jar
module purge
module use /apps/skl/modules/all
modList=("Java/17.0.6")
usage()
{
echo "# revertSam.sh Sort a BAM by read name then strip mapping info.  The result will be an unmapped BAM file.
# Dependencies:  Java v1.8+, GATK4
# Info: https://broadinstitute.github.io/picard/
#
# Usage: sbatch $0 -b /path/to/bam/folder -o /path/to/output/folder -S sampleID | [-h | --help]
#
# Options:
# -b <arg>           REQUIRED: Path to where your bam file is located
# -S <arg>           REQUIRED: ID of the sample which must be in the bam file name
# -o <arg>           Path to the output default: /hpcfs/users/$USER/sequences/sampleID/SLURM_JOB_ID
# -h | --help	     Prints the message you are reading.
#
# History:
# Script created by: Mark Corbett on 15/04/2020
# email:mark dot corbett is at adelaide university
# Modified (Date, Name, Description):
#
"
}

# Parse script options
while [ "$1" != "" ]; do
    case $1 in
        -b )            shift
                        bamDir=$1
                        ;;
        -S )            shift
                        sampleID=$1
                        ;;
        -o )            shift
                        outDir=$1
                        ;;
        -h | --help )   for mod in "${modList[@]}"; do
                            module load $mod
                        done
                        java -jar $PICARD RevertSam --help
                        for mod in "${modList[@]}"; do
                            module unload $mod
                        done
                        usage
                        exit 0
                        ;;
        * )             for mod in "${modList[@]}"; do
                            module load $mod
                        done
                        java -jar $PICARD RevertSam --help
                        for mod in "${modList[@]}"; do
                            module unload $mod
                        done
                        usage
                        exit 1
    esac
    shift
done

# Check that your script has everything it needs to start.
if [ -z "$bamDir" ]; then # If bamFile not specified then do not proceed
    usage
    echo "## ERROR: You need to specify -b /path/to/bamfile
# -b <arg>    REQUIRED: Path to where your bam file is located"
    exit 1
fi
if [ -z "$sampleID" ]; then # If sample not specified then do not proceed
    usage
    echo "## ERROR: You need to specify -S sampleID because I need this to make your file names
    # -S <arg>    ID of the sample which must be in the bam file name"
    exit 1
fi

bamFile=$( find $bamDir/*.bam | grep $sampleID )

if [ -z "$outDir" ]; then # If output directory not specified then make one up
    outDir=/hpcfs/users/${USER}/sequences/$sampleID/$SLURM_JOB_ID
    echo "## INFO: You didn't specify an output directory so I'm going to put your files here.
    $outDir"
fi
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi

# Load modules
for mod in "${modList[@]}"; do
    module load $mod
done

#Do the thing
java -Xmx16G -jar $PICARD RevertSam \
I=$bamFile \
O=$outDir/$sampleID.u.bam \
MAX_RECORDS_IN_RAM=10000000 \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=$outDir/tmp.$SLURM_JOB_ID
