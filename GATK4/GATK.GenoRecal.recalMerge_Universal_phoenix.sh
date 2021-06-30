#!/bin/bash

#SBATCH -J recal_merge
#SBATCH -o /hpcfs/users/%u/log/recalMerge-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --time=05:00:00
#SBATCH --mem=12GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

# Script that genotypes and refines variant calls on multiple samples
# Script variables (set and forget)
modList=("arch/haswell" "Java/1.8.0_121" "HTSlib/1.10.2-foss-2016b")

usage()
{
echo"
# Script that genotypes and refines variant calls on multiple samples
# Requires: GATKv4.x, Java, samtools
#
# Example: sbatch $0 -p prefix_for_files [-c /path/to/config -o /path/to/output] | [-h | --help]
#
# Options:
# -p Prefix               REQUIRED. Must specify the prefix used in the preceeding genDB script
# -o /path/to/output      OPTIONAL. Where you want to find your VCF and other files if not set then current directory will be used
# -c /path/to/config.cfg  OPTIONAL. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -h | --help             OPTIONAL. Displays this message
#
# Thomas Litster; 2/6/2021
#
# Modified: (Date; Name; Description)
#
"
}


## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -c )            shift
                                Config=$1
                                ;;
                -p )            shift
                                outPrefix=$1
                                ;;
                -o )            shift
                                workDir=$1
                                ;;
                -h | --help )   usage
                                exit 0
                                ;;
                * )             usage
                                exit 1
        esac
        shift
done
if [ -z "$Config" ]; then # If no Config file specified use the default
    Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config
if [ -z "$outPrefix" ]; then #If no outPrefix then fail
    usage
    echo "##ERROR: You need to specify the file prefix from the preceeding genDB script
# -p Prefix    REQUIRED. Must specify the prefix used in the preceeding genDB script"
    exit 1
fi

tmpDir=/hpcfs/groups/phoenix-hpc-neurogenetics/tmp/${USER}/$outPrefix # Use a tmp directory for all of the GATK temp files
if [ ! -d "$tmpDir" ]; then # If tmp directory doesn't exist then ask user to check prefix.
        echo "##ERROR: I can't locate $tmpDir. Is the path correct?"
        echo "The tmp directory is dependent on the name of the file prefix you specified check that this is correct."
        exit 1
fi
if [ -z "$workDir" ]; then #If workDir not specified then output to the default directory
        workDir=/hpcfs/groups/phoenix-hpc-neurogenetics/variants/vcf/$BUILD
        echo "Using $workDir as the output directory"
fi
if [ ! -d "$workDir" ]; then
    mkdir -p $workDir
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

cat $tmpDir/*.$outPrefix.pipeline.log >> $workDir/$outPrefix.pipeline.log
find $tmpDir/*.$outPrefix.all.recal.vcf > $tmpDir/$outPrefix.vcf.list.txt
sed 's,^,-I ,g' $tmpDir/$outPrefix.vcf.list.txt > $tmpDir/$outPrefix.input.vcf.list.txt

java -Xmx8g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar GatherVcfs  \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
$(cat ${tmpDir}/$outPrefix.input.vcf.list.txt) \
-O ${workDir}/${outPrefix}.${BUILD}.vcf >> ${workDir}/${outPrefix}.pipeline.log  2>&1

bgzip ${workDir}/${outPrefix}.${BUILD}.vcf
tabix ${workDir}/${outPrefix}.${BUILD}.vcf.gz

grep ERROR $workDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.pipeline.log > $workDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.pipeline.$outPrefix.ERROR.log

if [ -z $(cat $workDir/$outPrefix.pipeline.ERROR.log) ]; then
	rm $workDir/$outPrefix.pipeline.ERROR.log
    rm -r $tmpDir
else 
	echo "Some bad things went down while this script was running please see $outPrefix.pipeline.ERROR.log and prepare for disappointment."
	exit 1
fi
