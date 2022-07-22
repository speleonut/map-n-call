#!/bin/bash
#SBATCH -J liftOver
#SBATCH -o /hpcfs/users/%u/log/liftOver.slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch            	                            # partition (this is the queue your job will be added to)
#SBATCH -N 1                                                # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 4                                                # number of cores (here uses 4)
#SBATCH --time=08:00:00                                     # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=16G                                            # memory pool for all cores (here set to 16 GB)

# Notification configuration
#SBATCH --mail-type=END					    # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                		    # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  	    # Email to which notification will be sent

# liftOver.picard.sh
# Set location of picard.jar
PICARD=/hpcfs/groups/phoenix-hpc-neurogenetics/executables/gatk-latest/GenomeAnalysisTK.jar
modList=("arch/haswell" "Java/10.0.1")
usage()
{
echo "# liftOver.picard.sh Sort a BAM by read name then strip mapping info.  The result will be an unmapped BAM file.
# Dependencies:  Java v1.8+, GATK4
# Info: https://broadinstitute.github.io/picard/
#
# Usage: sbatch $0 -v /path/to/vcf/folder -C /path/to/chain -R /path/to/target_ref [ -B old_build,new_build -o /path/to/output/folder ] | [-h | --help]
# Usage: sbatch --array 0-(n-1 number of VCF) $0 -v /path/to/vcf/folder -C /path/to/chain -R /path/to/target_ref [ -o /path/to/output/folder ] | [-h | --help]
#
# Options:
# -v <arg>           REQUIRED: Path to where your vcf file is located
# -C <arg>           REQUIRED: Path to the chain file
# -R <arg>           REQUIRED: Path to reference sequence
# -o <arg>           OPTIONAL: Path to the output default: /hpcfs/users/$USER/liftOver/SLURM_JOB_ID
# -B <arg>,<arg>     OPTIONAL: If you have a genome build in your file name prior to .vcf e.g. *.hg19.vcf you can update to the new build e.g. *.hg38.vcf
# -h | --help	     Prints the message you are reading.
#
# History:
# Script created by: Mark Corbett on 22/07/2022 (3.143 <> pi day)
# email:mark dot corbett is at adelaide university
# Modified (Date, Name, Description):
#
"
}

# Parse script options
while [ "$1" != "" ]; do
    case $1 in
        -v )            shift
                        vcfDir=$1
                        ;;
        -C )            shift
                        chainFile=$1
                        ;;
        -R )            shift
                        targetGenome=$1
                        ;;
        -o )            shift
                        outDir=$1
                        ;;
        -B )            shift
                        oldBuild=${1%,*} # the bit before the comma
                        newBuild=${1##*,} # The bit after the comma
                        ;;
        -h | --help )   for mod in "${modList[@]}"; do
                            module load $mod
                        done
                        java -jar $PICARD LiftoverVcf --help
                        for mod in "${modList[@]}"; do
                            module unload $mod
                        done
                        usage
                        exit 0
                        ;;
        * )             for mod in "${modList[@]}"; do
                            module load $mod
                        done
                        java -jar $PICARD LiftoverVcf --help
                        for mod in "${modList[@]}"; do
                            module unload $mod
                        done
                        usage
                        exit 1
    esac
    shift
done

# Check that your script has everything it needs to start.
if [ -z "$vcfDir" ]; then # If vcf folder is not specified then do not proceed
    usage
    echo "## ERROR: You need to specify -v /path/to/vcf/folder
    # -v <arg>    REQUIRED: Path to where your vcf file or files are located"
    exit 1
fi
if [ -z "$chainFile" ]; then # If chain file is not specified then do not proceed
    usage
    echo "## ERROR: You need to specify -C /path/to/chain/file because I need this to do the lift over.  
    # You usally get the chain files from UCSC genome browser downloads
    # -C <arg>    REQUIRED: Path to the chain file"
    exit 1
fi
if [ -z "$targetGenome" ]; then # If genome file is not specified then do not proceed
    usage
    echo "## ERROR: You need to specify -R /path/to/target_genome because I need this to do the lift over.  
    # Check /hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq for available options
    # -R <arg>    REQUIRED: Path to the reference sequence"
    exit 1
fi
if [ -z "$oldBuild" ]; then
    vcfFile=($(find $vcfDir/*.vcf))
    vcfOut=$(basename ${vcfFile[$SLURM_ARRAY_TASK_ID]} .vcf).$(basename ${targetGenome}).vcf
else
    vcfFile=($(find $vcfDir/*.${oldBuild}.vcf))
    vcfOut=$(basename ${vcfFile[$SLURM_ARRAY_TASK_ID]} .${oldBuild}.vcf).${newBuild}.vcf
fi

if [ -z "$outDir" ]; then # If output directory not specified then make one up
    outDir=/hpcfs/users/${USER}/liftOver/$SLURM_JOB_ID
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
java -Xmx16G -jar picard.jar LiftoverVcf \
    I=${vcfFile[$SLURM_ARRAY_TASK_ID]} \
    O=$outDir/$vcfOut \
    CHAIN=$chainFile \
    REJECT=$outDir/rejected_variants.$vcfOut.gz \
    R=$targetGenome \
    CREATE_INDEX=true
