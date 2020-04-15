#!/bin/bash
# See usage for description and history

# Script variables (set and forget)
GATKPATH=/data/neurogenetics/executables/GenomeAnalysisTK-3.7 # Where the GATK program is even if in your $PATH
GATKREFPATH=/data/neurogenetics/RefSeq/GATK #Refseq index library location.
GATKINDEX=hg19_1stM_unmask_ran_all.fa
ChrIndexPath=$GATKREFPATH/$GATKINDEX.chridx #Location of index bed files
IndexBedFiles=01.hg19-M1.bed,02.hg19-2-3.bed,03.hg19-4-5.bed,04.hg19-6-7.bed,05.hg19-8-10.bed,06.hg19-11-13.bed,07.hg19-14-17.bed,08.hg19-18etc.bed # A comma separated array of names of index files
arrIndexBedFiles=$(echo $IndexBedFiles | tr "," "\n")
#SCRIPTPATH=/home/users/mcorbett/Scripts # Location of the cleanup script template
#BUILD=$(echo $GATKINDEX | awk '{print substr($1, 1, length($1) - 3)}') # Genome build used = $GATKINDEX less the .fa, this will be incorporated into file names.
BUILD=hg19 # If using hg19_1stM_unmask_ran_all.fa as GATKINDEX
tmpDir=$FASTDIR/tmp/$SLURM_JOBID # Use a tmp directory in $FASTDIR/tmp for all of the GATK temp files

usage()
{
echo"
# $0 [-p prefix_for_files -i /path/to/input_dir -o /path/to/output | -h | --help]
# Script that genotypes and refines variant calls on multiple samples then splits to individual VCF using GATK
# Requires: GATKv3.x, a folder of gVCF files, VCFSplitWithGATKv3.sh.
#
# Options:
# -p Prefix		Default is date(YYYYmmdd)-\$SLURM_JOBID. Prefix to help identify your merged gVCF files (This can be any alphanumeric text a short easily traced code is best eg. Project-Date)
# -i /path/to/input_dir	Where your gVCF are stored. Default is ~/DATA/gVcfFileLibrary
# -o /path/to/output	Where you want to find your VCF and other files if not set then current directory will be used
# -h | --help		Displays this message
#
# Example
# $0 -p MyProjectID -i ~/gVcfLibrary -o ~/Project-Output
#
# The current system variables set are:
# GATKPATH=$GATKPATH 		# Where the GATK program is
# GATKREFPATH=$GATKREFPATH 	# Refseq index library locations for GATK
# GATKINDEX=$GATKINDEX		# Base name of GATK indexes
# SCRIPTPATH=$SCRIPTPATH 	# Location of the cleanup script template
# IndexBedFiles=$IndexBedFiles
# BUILD=$BUILD				# Genome build code to be incorporated into file names
# tmpDir=$tmpDir			# Location of temp files will be deleted at the end of the script
#
# Forked from GATKv2.x.VCFRecalibration.HPC.pipeline.sh by Mark Corbett on 15/04/2014.
# mark.corbett@adelaide.edu.au
# Modified: (Date; Name; Description)
# 18/09/2014; Mark Corbett; Update to GATK 3.2-2
# 24/06/2015; Mark Corbett; Split initial genotyping step; Increase threads; Add bgzip; Auto OUTPREFIX; Remove option to put in an interval
# 08/10/2015; Mark Corbett; Update to GATKv3.4.x; Now search for vcf.gz$ files; Add toggle for Broad19.fasta index; Use different GATK bundle files based on BUILD
# 28/09/2016; Mark Corbett; Update to GATKv3.6
# 12/05/2017; Atma Ivancevic; Translated for use on phoenix HPC
"
}

# Set variables
while [ "$1" != "" ]; do
	case $1 in
		-p )			shift						
						OUTPREFIX=$1
						;;
		-i )			shift						
						inputDir=$1
						;;
		-o )			shift
						WORKDIR=$1
						;;
		-h | --help )			usage
						exit 1
						;;
		* )				usage
						exit 1
	esac
	shift
done
if [ -z "$OUTPREFIX" ]; then #If no outprefix specified then make one up
	OUTPREFIX=$(date "+%Y%m%d")"-"$SLURM_JOBID
	echo "Your files will be prefixed with the code: $OUTPREFIX"
fi
if [ -z "$inputDir" ]; then #If no input directory specified then use default
	inputDir=/data/neurogenetics/gVcfDumpingGround/Genomes/
	echo "Using $inputDir as the input directory"
fi
if [ -z "$WORKDIR" ]; then #If WORKDIR not specified then run in current directory
	WORKDIR=$(pwd)
	echo "Using current directory as the working directory"
fi
if [ ! -d $tmpDir ]; then # If tmp directory doesn't exist (likely) then create it.
	mkdir -p $tmpDir # -p forces the entire path to be created because if it doesn't exist epic failure ensues
fi

## Start script ##
mkdir -p $WORKDIR/$SLURM_JOBID
cd $inputDir

# Gather VCF to one file
find *.vcf.gz > $tmpDir/fileOfFiles.txt

cd $tmpDir
sed 's,^,--variant '"$inputDir"'\/,g' fileOfFiles.txt > inputGVCF.txt
rm fileOfFiles.txt

for bed in $arrIndexBedFiles; do
	mkdir -p $tmpDir/$bed
done

for bed in $arrIndexBedFiles; do
	java -Xmx20g -Djava.io.tmpdir=$tmpDir/$bed -jar $GATKPATH/GenomeAnalysisTK.jar \
	-T GenotypeGVCFs \
	-L $ChrIndexPath/$bed \
	-R $GATKREFPATH/$GATKINDEX \
	$(cat inputGVCF.txt) \
	-o $tmpDir/$bed.$OUTPREFIX.snps.vcf > $WORKDIR/$SLURM_JOBID/$bed.$OUTPREFIX.pipeline.log 2>&1 &
done
wait

cat $WORKDIR/$SLURM_JOBID/*.log >> $WORKDIR/$OUTPREFIX.pipeline.log
find *.$OUTPREFIX.snps.vcf > $OUTPREFIX.vcf.list.txt
sed 's,^,-V '"$tmpDir"'\/,g' $OUTPREFIX.vcf.list.txt > $OUTPREFIX.inputVCF.txt

java -cp $GATKPATH/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
-R $GATKREFPATH/$GATKINDEX \
-out $OUTPREFIX.snps.vcf \
$(cat $OUTPREFIX.inputVCF.txt) \
-assumeSorted >> $WORKDIR/$OUTPREFIX.pipeline.log  2>&1
  
# Variant quality recalibration
# Generate recalibration data for SNPs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R $GATKREFPATH/$GATKINDEX \
-input $OUTPREFIX.snps.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATKREFPATH/hapmap_3.3.$BUILD.vcf \
-resource:omni,known=false,training=true,truth=false,prior=12.0 $GATKREFPATH/1000G_omni2.5.$BUILD.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATKREFPATH/1000G_phase1.snps.high_confidence.$BUILD.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKREFPATH/dbsnp_138.$BUILD.vcf \
-an DP -an MQ -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP \
-nt 16 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile $OUTPREFIX.snp.recal \
-rscriptFile $OUTPREFIX.recal.snp.plots.R \
-tranchesFile $OUTPREFIX.snp.tranches >> $WORKDIR/$OUTPREFIX.pipeline.log 2>&1

# Apply recalibration for SNPs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R $GATKREFPATH/$GATKINDEX \
-input $OUTPREFIX.snps.vcf \
-mode SNP \
-recalFile $OUTPREFIX.snp.recal \
-tranchesFile $OUTPREFIX.snp.tranches \
--ts_filter_level 99.5 \
-nt 16 \
-o $OUTPREFIX.snps.recal.vcf >> $WORKDIR/$OUTPREFIX.pipeline.log 2>&1

#Generate recalibration data for INDELs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R $GATKREFPATH/$GATKINDEX \
-input $OUTPREFIX.snps.recal.vcf \
--maxGaussians 4 \
-resource:mills,known=true,training=true,truth=true,prior=12.0 $GATKREFPATH/Mills_and_1000G_gold_standard.indels.$BUILD.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKREFPATH/dbsnp_138.$BUILD.vcf \
-an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode INDEL \
-nt 16 \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-recalFile $OUTPREFIX.indel.recal \
-rscriptFile $OUTPREFIX.recal.indel.plots.R \
-tranchesFile $OUTPREFIX.indel.tranches >> $WORKDIR/$OUTPREFIX.pipeline.log 2>&1

# Apply recalibration for INDELs
java -Xmx30g -Djava.io.tmpdir=$tmpDir -jar $GATKPATH/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R $GATKREFPATH/$GATKINDEX \
-input $OUTPREFIX.snps.recal.vcf \
-mode INDEL \
-recalFile $OUTPREFIX.indel.recal \
-tranchesFile $OUTPREFIX.indel.tranches \
--ts_filter_level 99.0 \
-nt 16 \
-o $WORKDIR/$OUTPREFIX.all.recal.vcf >> $WORKDIR/$OUTPREFIX.pipeline.log 2>&1

cd $WORKDIR
bgzip $OUTPREFIX.all.recal.vcf
tabix $OUTPREFIX.all.recal.vcf.gz

grep ERROR $WORKDIR/$OUTPREFIX.pipeline.log > $WORKDIR/$OUTPREFIX.pipeline.ERROR.log
if [ -z $(cat $WORKDIR/$OUTPREFIX.pipeline.ERROR.log) ]; then
	rm -r $tmpDir
	rm $WORKDIR/$OUTPREFIX.pipeline.ERROR.log
else 
	echo "Some bad things went down while this script was running please see $OUTPREFIX.pipeline.ERROR.log and prepare for disappointment."
	rm -r $tmpDir
	exit 1
fi
