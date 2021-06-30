#!/bin/bash

#SBATCH -J VQSR
#SBATCH -o /hpcfs/users/%u/log/VQSR-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --time=24:00:00
#SBATCH --mem=36GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# See usage for description and history
# Script variables (set and forget)
modList=("arch/haswell" "Java/1.8.0_121" "arch/skylake" "R/4.0.3")

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

# Define files for the array
bedFile=($arrIndexBedFiles)

## Start script ##
cd ${tmpDir}
#Generate recalibration data for INDELs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar VariantRecalibrator \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.sites.only.vcf \
--max-gaussians 4 \
--resource:mills,known=true,training=true,truth=true,prior=12.0 $GATKREFPATH/$BUILD/$Mills_INDELS \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKREFPATH/$BUILD/$DBSNP \
--resource:axiomPoly,known=false,training=true,truth=false,prior=10 $GATKREFPATH/$BUILD/$Axiom \
-an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode INDEL \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-O ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.indel.recal \
--rscript-file ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.recal.indel.plots.R \
--tranches-file ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.indel.tranches >> $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.pipeline.log 2>&1

# Apply recalibration for INDELs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar ApplyVQSR \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.vcf \
-mode INDEL \
--recal-file ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.indel.recal \
--tranches-file ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.indel.tranches \
-ts-filter-level 99.0 \
-O $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.indel.recal.vcf >> $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.pipeline.log 2>&1

# Generate Recalibration data for SNPs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar VariantRecalibrator \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.sites.only.vcf \
--max-gaussians 6 \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATKREFPATH/$BUILD/$hapMap \
--resource:omni,known=false,training=true,truth=false,prior=12.0 $GATKREFPATH/$BUILD/$Omni \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATKREFPATH/$BUILD/${OneKg_HC_SNPs} \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKREFPATH/$BUILD/$DBSNP \
-an DP -an MQ -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.0 -tranche 90.0 \
-O ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.snp.recal \
--rscript-file ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.recal.snp.plots.R \
--tranches-file ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.snp.tranches >> $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.pipeline.log 2>&1

# Apply recalibration for SNPs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar ApplyVQSR \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.indel.recal.vcf \
-mode SNP \
--recal-file ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.snp.recal \
--tranches-file ${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.snp.tranches \
-ts-filter-level 99.5 \
-O $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.all.recal.vcf >> $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.pipeline.log 2>&1