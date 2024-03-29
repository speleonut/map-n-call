#!/bin/bash
# This is a coordinator script for genotyping of gVGFs with GATK4 and subsequent merging of the produced VCFs
whereAmI="$(dirname "$(readlink -f "$0")")" # Assumes that the script is linked to the git repo and the driectory structure is not broken
configDir="$(echo ${whereAmI} | sed -e 's,GATK4,configs,g')"
enviroCfg="${configDir}/BWA-GATKHC.environment.cfg"
source ${enviroCfg}

if [ ! -d "${logDir}" ]; then
    mkdir -p ${logDir}
    echo "## INFO: New log directory created, you'll find all of the log information from this pipeline here: ${logDir}"
fi

usage()
{
echo "
# Coordinator script for genotype recalibration and subsequent VCF merging
#
# Options:
# -p Prefix                     OPTIONAL. Default is date(YYYYmmdd_unix_time). Prefix to help identify your merged gVCF files (This can be any alphanumeric text a short easily traced code is best eg. Project-Date)
# -i /path/to/sample-name-map   RECOMMENDED. Location of a list of sample names and gVCF locations in a tab delimited format.  If not supplied the script will scrape the gVcfFolder as specified in the Config file.
# -o /path/to/output            OPTIONAL. Where you want to find your VCF and other files if not set then current directory will be used
# -c /path/to/Config.cfg        OPTIONAL. A default Config will be used if this is not specified.  The Config contains all of the stuff that used to be set in the top part of our scripts
# --gpfs                        OPTIONAL. Use /gpfs/users/$USER/tmp as your tmp directory for better read/write performance.
# -h | --help                   OPTIONAL. Displays this message
#
# Example:
# $0 -p MyProjectID -i ~/gVcfLibrary -o ~/Project-Output
#
# Thomas Litster; 2/6/2021
#
# Modified: (Date; Name; Description)
# 07/07/2021; Mark Corbett; Update to accept any genome build, change default output paths
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
                --gpfs )        shift
                                tmpDir="/gpfs/users/${USER}/tmp"
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

if [ -z "$outPrefix" ]; then #If no outPrefix specified then make one up
    outPrefix=$(date "+%Y%m%d_%s")
    echo "## INFO: Your VCF files will be prefixed with the code: $outPrefix"
fi

tmpDir=$tmpDir/$outPrefix # Use a tmp directory for all of the GATK and samtools temp files
echo "## INFO: Using ${tmpDir} as the tmp directory location."
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi
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

## Submit Jobs ##
genotypeJob=`sbatch --array=0-23 --export=ALL,enviroCfg=${enviroCfg},tmpDir=${tmpDir} $scriptDir/GATK4/GATK.GenoRecal.GenDBGeno_Universal_phoenix.sh -c ${Config} -i ${sampleNameMap} -o ${workDir} -p ${outPrefix}`
genotypeJob=$(echo $genotypeJob | cut -d" " -f4)
sbatch --export=ALL,enviroCfg=${enviroCfg},tmpDir=${tmpDir} --dependency=afterok:${genotypeJob} $scriptDir/GATK4/GATK.GenoRecal.recalMerge_Universal_phoenix.sh -c ${Config} -p ${outPrefix} -o ${workDir}
