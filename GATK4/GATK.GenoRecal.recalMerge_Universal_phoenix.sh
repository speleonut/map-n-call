#!/bin/bash

#SBATCH -J recal_merge
#SBATCH -o /hpcfs/users/%u/log/recalMerge-slurm-%j.out
#SBATCH -p skylake,icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --time=06:00:00
#SBATCH --mem=48GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

# Script that genotypes and refines variant calls on multiple samples
# Script variables (set and forget)
source ${enviroCfg}
module purge
module use /apps/skl/modules/all
modList=("R/4.3.1-foss-2021b" "Java/17.0.6" "HTSlib/1.17-GCC-11.2.0")

usage()
{
echo "
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

if [ -z "$outPrefix" ]; then #If no outPrefix then fail
    usage
    echo "##ERROR: You need to specify the file prefix from the preceeding genDB script"
    echo "# -p Prefix    REQUIRED. Must specify the prefix used in the preceeding genDB script"
    exit 1
fi

if [ -z "${scriptDir}" ]; then # Test if the script was executed independently of the Universal Launcher script
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,GATK4,configs,g')"
    source ${configDir}/BWA-GATKHC.environment.cfg
    if [ ! -d "${logDir}" ]; then
        mkdir -p ${logDir}
        echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
    fi
    tmpDir=${tmpDir}/${outPrefix}
    if [ ! -d "$tmpDir" ]; then # If tmp directory doesn't exist then ask user to check prefix.
        echo "##ERROR: I can't locate $tmpDir. Is the path correct?"
        echo "The tmp directory is dependent on the name of the file prefix you specified and it should exist or the is no point in running this script, check that this is correct."
        exit 1
    fi
fi

if [ -z "$Config" ]; then # If no Config file specified use the default
    Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config

if [ -z "$workDir" ]; then #If workDir not specified then output to the default directory
        workDir=/hpcfs/groups/phoenix-hpc-neurogenetics/variants/vcf/$BUILD/$outPrefix
        echo "Using $workDir as the output directory"
fi
if [ ! -d "$workDir" ]; then
    mkdir -p $workDir
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

#Collect all log files
cat $tmpDir/*.${outPrefix}.${BUILD}.pipeline.log >> $workDir/${outPrefix}.${BUILD}.pipeline.log

#Merge and sort VCFs from the previous step
find $tmpDir/*.$outPrefix.sites.only.vcf > $tmpDir/$outPrefix.vcf.list.txt
sed 's,^,-I ,g' $tmpDir/$outPrefix.vcf.list.txt > $tmpDir/$outPrefix.input.vcf.list.txt
$GATKPATH/gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmpDir" SortVcf  \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
--MAX_RECORDS_IN_RAM 2000000 \
--TMP_DIR $tmpDir \
$(cat ${tmpDir}/$outPrefix.input.vcf.list.txt) \
-O ${tmpDir}/${outPrefix}.merge.sites.only.vcf >> ${workDir}/${outPrefix}.${BUILD}.pipeline.log  2>&1

find $tmpDir/*.$outPrefix.vcf > $tmpDir/$outPrefix.vcf.list.txt
sed 's,^,-I ,g' $tmpDir/$outPrefix.vcf.list.txt > $tmpDir/$outPrefix.input.vcf.list.txt
$GATKPATH/gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmpDir" SortVcf  \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
--MAX_RECORDS_IN_RAM 2000000 \
--TMP_DIR $tmpDir \
$(cat ${tmpDir}/$outPrefix.input.vcf.list.txt) \
-O ${tmpDir}/${outPrefix}.merge.vcf >> ${workDir}/${outPrefix}.${BUILD}.pipeline.log  2>&1

#Generate recalibration data for INDELs
$GATKPATH/gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmpDir" VariantRecalibrator \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/${outPrefix}.merge.sites.only.vcf \
--max-gaussians 4 \
--resource:mills,known=true,training=true,truth=true,prior=12.0 $GATKREFPATH/$BUILD/$Mills_INDELS \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKREFPATH/$BUILD/$DBSNP \
--resource:axiomPoly,known=false,training=true,truth=false,prior=10 $GATKREFPATH/$BUILD/$Axiom \
-an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode INDEL \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-O $tmpDir/$outPrefix.indel.recal \
--rscript-file $workDir/$outPrefix.recal.indel.plots.R \
--tranches-file $tmpDir/$outPrefix.indel.tranches >> $workDir/${outPrefix}.${BUILD}.pipeline.log 2>&1

# Apply recalibration for INDELs
$GATKPATH/gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmpDir" ApplyVQSR \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/$outPrefix.merge.vcf \
-mode INDEL \
--recal-file $tmpDir/$outPrefix.indel.recal \
--tranches-file $tmpDir/$outPrefix.indel.tranches \
-ts-filter-level 99.0 \
-O $tmpDir/$outPrefix.indel.recal.vcf >> $workDir/${outPrefix}.${BUILD}.pipeline.log 2>&1

# Generate Recalibration data for SNPs
$GATKPATH/gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmpDir" VariantRecalibrator \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/${outPrefix}.merge.sites.only.vcf \
--max-gaussians 4 \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATKREFPATH/$BUILD/$hapMap \
--resource:omni,known=false,training=true,truth=false,prior=12.0 $GATKREFPATH/$BUILD/$Omni \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATKREFPATH/$BUILD/${OneKg_HC_SNPs} \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKREFPATH/$BUILD/$DBSNP \
-an DP -an MQ -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.0 -tranche 90.0 \
-O $tmpDir/$outPrefix.snp.recal \
--rscript-file $workDir/$outPrefix.recal.snp.plots.R \
--tranches-file $tmpDir/$outPrefix.snp.tranches >> $workDir/${outPrefix}.${BUILD}.pipeline.log 2>&1

# Apply recalibration for SNPs
$GATKPATH/gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmpDir" ApplyVQSR \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/$outPrefix.indel.recal.vcf \
-mode SNP \
--recal-file $tmpDir/$outPrefix.snp.recal \
--tranches-file $tmpDir/$outPrefix.snp.tranches \
-ts-filter-level 99.5 \
-O $workDir/$outPrefix.${BUILD}.vcf >> $workDir/${outPrefix}.${BUILD}.pipeline.log 2>&1

bgzip ${workDir}/${outPrefix}.${BUILD}.vcf
tabix ${workDir}/${outPrefix}.${BUILD}.vcf.gz

## Check for bad things and clean up
if [ ! -f "${workDir}/${outPrefix}.${BUILD}.vcf.gz.tbi" ]; then
    echo "##ERROR: Some bad things went down while this script was running please see $workDir/${outPrefix}.${BUILD}.pipeline.ERROR.log and prepare for disappointment."
    exit 1
fi

grep -i ERROR $workDir/${outPrefix}.${BUILD}.pipeline.log > $workDir/${outPrefix}.${BUILD}.pipeline.ERROR.log
if [ -z "$(cat $workDir/${outPrefix}.${BUILD}.pipeline.ERROR.log)" ]; then
	rm $workDir/${outPrefix}.${BUILD}.pipeline.ERROR.log
    rm -r $tmpDir
else 
	echo "Some bad things went down while this script was running please see $workDir/${outPrefix}.${BUILD}.pipeline.ERROR.log and prepare for disappointment."
	exit 1
fi
