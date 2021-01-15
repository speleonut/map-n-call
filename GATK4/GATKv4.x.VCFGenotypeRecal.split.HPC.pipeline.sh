#!/bin/bash

#SBATCH -J GATKGeno
#SBATCH -o /hpcfs/users/%u/launch/GATKgeno-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=12:00:00
#SBATCH --mem=64GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# See usage for description and history
# Script variables (set and forget)
tmpDir=/hpcfs/users/${USER}/tmp/$SLURM_JOBID # Use a tmp directory in $FASTDIR/tmp for all of the GATK temp files

usage()
{
echo"
# $0 [-p prefix_for_files -i /path/to/input_dir -o /path/to/output | -h | --help]
# Script that genotypes and refines variant calls on multiple samples then splits to individual VCF using GATK
# Requires: GATKv4.x, Java, samtools
#
# Options:
# -p Prefix		          OPTIONAL. Default is date(YYYYmmdd)-\$SLURM_JOBID. Prefix to help identify your merged gVCF files (This can be any alphanumeric text a short easily traced code is best eg. Project-Date)
# -i /path/to/input_dir	  OPTIONAL. Where your gVCF are stored. Default is $FASTDIR/VCF/gVCF/\$BUILD/genomes
# -o /path/to/output	  OPTIONAL. Where you want to find your VCF and other files if not set then current directory will be used
# -c /path/to/config.cfg. OPTIONAL. A default config will be used if this is not specified.  The config contains all of the stuff that used to be set in the top part of our scripts
# -h | --help		      OPTIONAL. Displays this message
#
# Example
# $0 -p MyProjectID -i ~/gVcfLibrary -o ~/Project-Output
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
#
"
}

# Set variables
while [ "$1" != "" ]; do
	case $1 in
		-c )			shift
					    config=$1
					    ;;
		-p )			shift						
						outPrefix=$1
						;;
		-i )			shift						
						inputDir=$1
						;;
		-o )			shift
						workDir=$1
						;;
		-h | --help )			usage
						exit 1
						;;
		* )				usage
						exit 1
	esac
	shift
done
if [ -z "$config" ]; then # If no config file specified use the default
   config=/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/GATK4/BWA-GATKHC.GRCh38_full_analysis_set_phoenix.cfg
fi
source $config

if [ -z "$outPrefix" ]; then #If no outPrefix specified then make one up
	outPrefix=$(date "+%Y%m%d")"-"$SLURM_JOBID
	echo "Your files will be prefixed with the code: $outPrefix"
fi
if [ -z "$inputDir" ]; then #If no input directory specified then use default
	inputDir=/hpcfs/users/${USER}/gVCF/
	echo "Using $inputDir as the input directory"
fi
if [ -z "$workDir" ]; then #If workDir not specified then run in current directory
	workDir=$(pwd)
	echo "Using current directory as the working directory"
fi
if [ ! -d "$tmpDir" ]; then # If tmp directory doesn't exist (likely) then create it.
	mkdir -p $tmpDir # -p forces the entire path to be created because if it doesn't exist epic failure ensues
fi

# load modules
module load arch/haswell
module load Java/1.8.0_121
module load HTSlib/1.3.1-GCC-5.3.0-binutils-2.25
module load R/3.4.2-foss-2016b
module load GCC/5.4.0-2.26 

## Start script ##
mkdir -p $workDir/$SLURM_JOBID
cd $inputDir

# Gather VCF to one file
find *.vcf.gz > $tmpDir/fileOfFiles.txt

cd $tmpDir
sed 's,^, -V '"$inputDir"'\/,g' fileOfFiles.txt > inputGVCF.txt
#rm fileOfFiles.txt

for bed in $arrIndexBedFiles; do
	mkdir -p $tmpDir/$bed
done

java -Xmx20g -Djava.io.tmpdir=$tmpDir/$bed -jar $GATKPATH/GenomeAnalysisTK.ja GenomicsDBImport \
$(cat inputGVCF.txt) \
--genomicsdb-workspace-path my_database \
--intervals $ChrIndexPath/$bed \


for bed in $arrIndexBedFiles; do
	java -Xmx20g -Djava.io.tmpdir=$tmpDir/$bed -jar $GATKPATH/GenomeAnalysisTK.jar GenotypeGVCFs \
	-L $ChrIndexPath/$bed \
	-R $GATKREFPATH/$BUILD/$GATKINDEX \
	-V $inputDir/FR07959028.snps.g.vcf.gz \
	-O $tmpDir/$bed.$outPrefix.snps.vcf > $workDir/$SLURM_JOBID/$bed.$outPrefix.pipeline.log 2>&1 &
done
wait

cat $workDir/$SLURM_JOBID/*.log >> $workDir/$outPrefix.pipeline.log
find *.$outPrefix.snps.vcf > $outPrefix.vcf.list.txt
sed 's,^,-I '"$tmpDir"'\/,g' $outPrefix.vcf.list.txt > $outPrefix.inputVCF.txt

java -Xmx8g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar GatherVcfs \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
$(cat $outPrefix.inputVCF.txt) \
-O $outPrefix.snps.vcf > $workDir/$outPrefix.pipeline.log  2>&1
  
# Variant quality recalibration
# Generate recalibration data for SNPs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar VariantRecalibrator \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/$outPrefix.snps.vcf \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATKREFPATH/$BUILD/hapmap_3.3.hg38.vcf \
--resource:omni,known=false,training=true,truth=false,prior=12.0 $GATKREFPATH/$BUILD/1000G_omni2.5.hg38.vcf \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATKREFPATH/$BUILD/1000G_phase1.snps.high_confidence.hg38.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKREFPATH/$BUILD/dbsnp_146.hg38.vcf \
-an DP -an MQ -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-O $outPrefix.snp.recal \
--rscript-file $outPrefix.recal.snp.plots.R \
--tranches-file $outPrefix.snp.tranches >> $workDir/$outPrefix.pipeline.log 2>&1

# Apply recalibration for SNPs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar ApplyVQSR \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/$outPrefix.snps.vcf \
-mode SNP \
--recal-file $outPrefix.snp.recal \
--tranches-file $outPrefix.snp.tranches \
-ts-filter-level 99.5 \
-O $outPrefix.snps.recal.vcf >> $workDir/$outPrefix.pipeline.log 2>&1

#Generate recalibration data for INDELs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar VariantRecalibrator \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/$outPrefix.snps.recal.vcf \
--max-gaussians 4 \
--resource:mills,known=true,training=true,truth=true,prior=12.0 $GATKREFPATH/$BUILD/Mills_and_1000G_gold_standard.indels.hg38.vcf \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKREFPATH/$BUILD/dbsnp_146.hg38.vcf \
--resource:axiomPoly,known=false,training=true,truth=false,prior=10 $GATKREFPATH/$BUILD/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf \
-an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode INDEL \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-O $outPrefix.indel.recal \
--rscript-file $outPrefix.recal.indel.plots.R \
--tranches-file $outPrefix.indel.tranches >> $workDir/$outPrefix.pipeline.log 2>&1

# Apply recalibration for INDELs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar ApplyVQSR \
-R $GATKREFPATH/$BUILD/$GATKINDEX \
-V $tmpDir/$outPrefix.snps.recal.vcf \
-mode INDEL \
--recal-file $outPrefix.indel.recal \
--tranches-file $outPrefix.indel.tranches \
-ts-filter-level 99.0 \
-O $workDir/$outPrefix.all.recal.vcf >> $workDir/$outPrefix.pipeline.log 2>&1

cd $workDir
bgzip $outPrefix.all.recal.vcf
tabix $outPrefix.all.recal.vcf.gz

grep ERROR $workDir/$outPrefix.pipeline.log > $workDir/$outPrefix.pipeline.ERROR.log
if [ -z $(cat $workDir/$outPrefix.pipeline.ERROR.log) ]; then
	rm -r $tmpDir
	rm $workDir/$outPrefix.pipeline.ERROR.log
else 
	echo "Some bad things went down while this script was running please see $outPrefix.pipeline.ERROR.log and prepare for disappointment."
	rm -r $tmpDir
	exit 1
fi
