#!/bin/bash

#SBATCH -J GATKGeno
#SBATCH -o /hpcfs/users/%u/log/genDBGeno-slurm-%j.out
#SBATCH -p skylake,icelake,a100cpu
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
source ${enviroCfg}
module purge
module use /apps/skl/modules/all
modList=("Java/17.0.6")

usage()
{
echo "
# Script that genotypes and refines variant calls on multiple samples then splits to individual VCF using GATK
# Requires: GATKv4.x, Java, samtools
#
# Example:
# sbatch $0 [-p prefix_for_files -i /path/to/sample-name-map -o /path/to/output -c /path/to/config.cg] | [-h | --help]
#
# Options:
# -p Prefix                     OPTIONAL. Default is date(YYYYmmdd_unix_time). Prefix to help identify your merged gVCF files (This can be any alphanumeric text a short easily traced code is best eg. Project-Date)
# -i /path/to/sample-name-map   RECOMMENDED. Location of a list of sample names and gVCF locations in a tab delimited format.  If not supplied the script will scrape the gVcfFolder as specified in the Config file.
# -o /path/to/output            OPTIONAL. Where you want to find your VCF and other files if not set then current directory will be used
# -c /path/to/Config.cfg        OPTIONAL. A default Config will be used if this is not specified.  The Config contains all of the stuff that used to be set in the top part of our scripts
# -h | --help                   OPTIONAL. Displays this message
#
# Forked from GATKv2.x.VCFRecalibration.HPC.pipeline.sh by Mark Corbett on 15/04/2014.
# mark.corbett@adelaide.edu.au
# Modified: (Date; Name; Description)
# 18/09/2014; Mark Corbett; Update to GATK 3.2-2
# 24/06/2015; Mark Corbett; Split initial genotyping step; Increase threads; Add bgzip; Auto outPrefix; Remove option to put in an interval
# 08/10/2015; Mark Corbett; Update to GATKv3.4.x; Now search for vcf.gz$ files; Add toggle for Broad19.fasta index; Use different GATK bundle files based on BUILD
# 28/09/2016; Mark Corbett; Update to GATKv3.6
# 12/05/2017; Atma Ivancevic; Translated for use on phoenix HPC
# 16/04/2020; Mark Corbett; Incorporate queue script GATKGenoRecal.split.q, incorporate config files.
# 28/2/2021; Mark Corbett,Ali Gardner; incorporate GenomicDB
# 2/5/2021; Thomas Litster; Corrected array submission and added vcf list to allow for subsequent merging
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
                -i )            shift
                                sampleNameMap=$1
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

if [ -z "$outPrefix" ]; then #If no outPrefix specified then make one up
    outPrefix=$(date "+%Y%m%d_%s")
    echo "## INFO: Your VCF files will be prefixed with the code: $outPrefix"
fi

if [ -z "${enviroCfg}" ]; then # Test if the script was executed independently of the Universal Launcher script
    whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
    configDir="$(echo ${whereAmI} | sed -e 's,GATK4,configs,g')"
    source ${configDir}/BWA-GATKHC.environment.cfg
    if [ ! -d "${logDir}" ]; then
        mkdir -p ${logDir}
        echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
    fi
    tmpDir=${tmpDir}/${outPrefix}
    if [ ! -d "$tmpDir" ]; then
        mkdir -p $tmpDir
    fi
fi

if [ -z "$Config" ]; then # If no Config file specified use the default
    Config=$scriptDir/configs/BWA-GATKHC.hs38DH_phoenix.cfg
    echo "## INFO: Using the default config ${Config}"
fi
source $Config

if [ -z "$workDir" ]; then #If workDir not specified then output to the default directory
        workDir=/hpcfs/groups/phoenix-hpc-neurogenetics/variants/vcf/$BUILD
        echo "## INFO: Using $workDir as the output directory"
fi
if [ ! -d "$workDir" ]; then
    mkdir -p $workDir
fi
if [ -z "$sampleNameMap" ]; then #If sampleNameMap is not supplied try to find some gVCFs to genotype
    echo "## WARN: A list of samples was not supplied so I'm going to genotype everything in $gVcfFolder."
    find $gVcfFolder/*.g.vcf.gz | cut -f1 -d"." | awk -F"/" '{print $NF}' > $tmpDir/$outPrefix.sample.list
    find $gVcfFolder/*.g.vcf.gz > $tmpDir/$outPrefix.vcf.list
    paste $tmpDir/$outPrefix.sample.list $tmpDir/$outPrefix.vcf.list > $workDir/$outPrefix.sample.name.map.txt
    sampleNameMap=$workDir/$outPrefix.sample.name.map.txt
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

# Define files for the array
bedFile=($arrIndexBedFiles)

# If this is a rerun GATK4 won't overwrite the previous database and will produce an error in the log
# therefore if there is a previous version it's best to remove it
if [ -d "$tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.${outPrefix}\_database" ]; then
    rm -rf $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.${outPrefix}\_database
    echo "## WARN: Possible script re-run detected, a previous genomics database $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.${outPrefix}\_database has been removed."
fi

## Start script ##
cd ${tmpDir}

$GATKPATH/gatk --java-options "-Xmx32g -Xms32g -Djava.io.tmpdir=$tmpDir" GenomicsDBImport \
--genomicsdb-workspace-path $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.${outPrefix}\_database \
--genomicsdb-shared-posixfs-optimizations true \
--batch-size 50 \
--reader-threads 4 \
--merge-input-intervals \
--consolidate \
--sample-name-map ${sampleNameMap} \
--intervals $ChrIndexPath/${bedFile[$SLURM_ARRAY_TASK_ID]} > $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.${outPrefix}.${BUILD}.pipeline.log  2>&1

$GATKPATH/gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmpDir" GenotypeGVCFs \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-D $GATKREFPATH/$BUILD/$DBSNP \
-G AS_StandardAnnotation \
-V gendb://${bedFile[$SLURM_ARRAY_TASK_ID]}.${outPrefix}\_database \
-O $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.vcf \
--merge-input-intervals >> $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.${outPrefix}.${BUILD}.pipeline.log  2>&1

$GATKPATH/gatk --java-options "-Xmx32g -Djava.io.tmpdir=$tmpDir" MakeSitesOnlyVcf \
-I $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.vcf \
-O $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.$outPrefix.sites.only.vcf >> $tmpDir/${bedFile[$SLURM_ARRAY_TASK_ID]}.${outPrefix}.${BUILD}.pipeline.log  2>&1
