#!/bin/bash
#SBATCH -J MAS
#SBATCH -o /fast/users/%u/log/revertSam.slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch            	                            # partition (this is the queue your job will be added to)
#SBATCH -N 1                                                # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 2                                                # number of cores (here uses 2)
#SBATCH --time=08:00:00                                     # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=8G                                            # memory pool for all cores (here set to 8 GB)

# Notification configuration
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                		    # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	    # Email to which notification will be sent

# revertSam.sh
usage()
{
echo "# revertSam.sh Sort a BAM by read name then strip mapping info.  The result will be an unmapped BAM file.
# Dependencies:  Java v1.8+, Picard v2.22.3+
# Info: https://broadinstitute.github.io/picard/
#
# Usage: sbatch $0 -b /path/to/bam/folder -o /path/to/output/folder -S sampleID | [-h | --help]
#
# Options:
# -b <arg>           REQUIRED: Path to where your bam file is located
# -S <arg>           REQUIRED: ID of the sample which must be in the bam file name
# -o <arg>           Path to the output default: $FASTDIR/sequences/sampleID/SLURM_JOB_ID
# -h | --help	     Prints the message you are reading.
#
# History:
# Script created by: Mark Corbett on 15/04/2020
# email:mark dot corbett is at adelaide university
# Modified (Date, Name, Description):
#
"
}
# Set location of picard.jar
PICARD=/data/neurogenetics/executables/Picard-latest/picard.jar

# Parse script options
while [ "$1" != "" ]; do
	case $1 in
		 -b )	shift
			bamDir=$1
 			;;
		 -S )	shift
			sampleID=$1
 			;;
		 -o )	shift
			outDir=$1
 			;;
		 -h | --help )	module load Java/10.0.1
		                java -jar $PICARD RevertSam --help
		                module unload Java/10.0.1
				            usage
				            exit 0
     ;;
		* )	module load Java/10.0.1
		    java -jar $PICARD RevertSam --help
			  module unload Java/10.0.1
			  usage
			  exit 1
	esac
	shift
done

# Check that your script has everything it needs to start.
if [ -z "$bamDir" ]; then # If bamFile not specified then do not proceed
    usage
    echo "#ERROR: You need to specify -b /path/to/bamfile
    # -b <arg>    REQUIRED: Path to where your bam file is located"
    exit 1
fi
if [ -z "$sampleID" ]; then # If sample not specified then do not proceed
    usage
    echo "#ERROR: You need to specify -S sampleID because I need this to make your file names
    # -S <arg>    ID of the sample which must be in the bam file name"
    exit 1
fi

bamFile=$( find $bamDir/*.bam | grep $sampleID )

if [ -z "$outDir" ]; then # If output directory not specified then make one up
    outDir=$FASTDIR/sequences/$sampleID/$SLURM_JOB_ID
    echo "#INFO: You didn't specify an output directory so I'm going to put your files here.
    $outDir"
fi
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi

# Load modules
module load Java/10.0.1

#Do the thing
java -Xmx4G -jar $PICARD \
-I=$bamFile \
-O=$outDir/$sampleID.u.bam \
--MAX_RECORDS_IN_RAM=4000000 \
--VALIDATION_STRINGENCY=LENIENT
