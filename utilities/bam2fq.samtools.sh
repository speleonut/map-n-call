#!/bin/bash
#SBATCH -J BAM2fq
#SBATCH -o /hpcfs/users/%u/log/bam2fq.samtools.slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1                            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 12                            # number of cores (here uses 12)
#SBATCH --time=05:30:00                 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=48G                       # memory pool for all cores (here set to 32 GB)

# Notification configuration
#SBATCH --mail-type=END                 # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  # Email to which notification will be sent

#Script Paths
userDir=/hpcfs/users/${USER}
refDir="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq"

# Script settings
modList=("SAMtools/1.17-GCC-11.2.0")
nCores=10

# Script functions
select_genome_build()
{
case "${genomeSize}" in
    3099922541 )    buildID="GRCh38"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
                    ;;
    3217346917 )    buildID="hs38DH"
                    genomeBuild="$refDir/hs38DH.fa"
                    ;;
    3137454505 )    buildID="hs37d5"
                    genomeBuild="$refDir/hs37d5.fa.gz"
                    ;;
    2730871774 )    buildID="GRCm38"   
                    genomeBuild="$refDir/GRCm38_68.fa"
                    ;;
    3117463893 )    buildID="CHM13v2"
                    genomeBuild="$refDir/T2T_CHM13v2.0.ucsc.ebv.fa.gz"
                    ;;
    3137161264 )    buildID="hg19"
                    genomeBuild="$refDir/ucsc.hg19.fasta"
                    ;;
    3105715063 )    buildID="GRCh38.hs38d1.no_alt"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
                    ;;
    3215250450 )    buildID="GRCh38.hs38d1.full"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz"
                    ;;
    3099750718 )    buildID="GRCh38"
                    genomeBuild="$refDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
                    ;;
    3031042417 )    buildID="GRCh38.blacklist"
                    genomeBuild="$refDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz"
                    ;;
    3101804741 )    buildID="hg19_1stM_unmask_ran_all"
                    genomeBuild="$refDir/hg19_1stM_unmask_ran_all.fa"
                    ;;
    * )         echo "## ERROR: Genome length $genomeSize for ${bamFile[SLURM_ARRAY_TASK_ID]} was not matched, you may need to specify the genome build directly using the -g flag."
                exit 1
                ;;
esac
echo "## INFO: The CRAM file ${bamFile[SLURM_ARRAY_TASK_ID]} was likely mapped to $buildID corresponding to the refseq $genomeBuild."
}

usage()
{
echo "# bam2fq.samtools.sh Sort a BAM by read name then convert to gzipped fastq files.
# Dependencies:  samtools v1.9+
# Info: http://www.htslib.org/doc/samtools.html
#
# Usage: sbatch $0 -b /path/to/bam/file.bam -o /path/to/output/folder -S sampleID | [-h | --help]
#        sbatch --array <0-(n-1) samples> $0 -i /path/to/input-file -o /path/to/output/folder
#
# Options:
# -b <arg>           REQUIRED: Path to your bam file
# -S <arg>           OPTIONAL: ID of the sample which will form the first part of your fastq file names. If not specified the sample name from the BAM file header will be used.
# -o <arg>           OPTIONAL: Path to the output. Note: the script will add the sample ID to this. Default: $userDir/sequences/bam2fq/sampleID
# -h | --help        Prints the message you are reading.
#
# Array job options:
# -i <arg>           REQUIRED: If running an array job a tab delimited text file with two columns listing the required arguments -b and -S from above (in that order).
# -o <arg>           OPTIONAL: Path to the output. Note: the script will add the sample ID to this. default: $userDir/sequences/bam2fq/sampleID
# -h | --help        Prints the message you are reading.
#
# History:
# Script created by: Mark Corbett on 15/04/2020
# email:mark dot corbett is at adelaide university
# Modified (Date; Name; Description):
# 15/01/2021; Mark Corbett; Update for new Phoenix and create array job option
# 07/12/2023; Mark Corbett; Change behaviour of -b to now specifiy the bam file. Increase memory allocation
#
"
}

# Parse script options
while [ "$1" != "" ]; do
    case $1 in
        -b )    shift
                bamFile=$1
                ;;
        -S )    shift
                sampleID=$1
                ;;
        -o )    shift
                outDir=$1
                ;;
        -i )    shift
                inputFile=$1
                ;;
        -h | --help )   for mod in "${modList[@]}"; do
                            module load $mod
                        done
                        samtools fastq
                        for mod in "${modList[@]}"; do
                            module unload $mod
                        done
                        usage
                        exit 0
                        ;;
        * )     for mod in "${modList[@]}"; do
                    module load $mod
                done
                samtools fastq
                for mod in "${modList[@]}"; do
                    module unload $mod
                done
                usage
                exit 1
    esac
    shift
done

# Check that your script has everything it needs to start.
if [ ! -z "$inputFile" ]; then
    bamFile=($(awk '{print $1}' $inputFile))
    sampleID=($(awk '{print $2}' $inputFile))
fi	
if [ -z "$bamFile" ]; then # If bamFile not specified then do not proceed
    usage
    echo "## ERROR: You need to specify -b /path/to/bam/file.bam
    # -b <arg>    REQUIRED: Path to where your bam file is located"
    exit 1
fi
# Fetch the genome build by looking at the bam header and then use that to select the correct reference sequence.
genomeSize=$(samtools view -H ${bamFile[SLURM_ARRAY_TASK_ID]} | grep @SQ | cut -f3 | cut -f2 -d":" | awk '{s+=$1} END {printf "%.0f\n", s}' -)
select_genome_build

# Check the sample ID
if [ -z "$sampleID" ]; then # If sample not specified then extract from the bam header
    sampleID=$(samtools view -H ${bamFile[SLURM_ARRAY_TASK_ID]} | grep "^@RG" | sed 's/.*SM:\([^[:space:]]\+\).*/\1/' | head -n 1)
    echo "## INFO: No sample ID provided. Sample ID ${sampleID} extracted from the bam header will be used to name the output files."
    if [ -z "$sampleID" ]; then # You've got problems!
        echo "## ERROR: No sample ID provided and could not extract a sample ID from the bam header. Please check your bam file and/or provide a sample ID using the -S option."
        usage
        exit 1
    fi
else
    sampleID=${sampleID[SLURM_ARRAY_TASK_ID]} # Capture the sample ID for this array job from the input file (or if was specified as an option then this won't change it).
fi

# Check the output directory and if it wasn't provided create a default directory.
if [ -z "${outDir}" ]; then
    outDir=$userDir/sequences/bam2fq/${sampleID}
    echo "## INFO: You didn't specify an output directory so I'm going to put your files here.
    BAM2FQ_DIR:${outDir}"
fi
if [ ! -d "${outDir}" ]; then
    mkdir -p ${outDir}/${sampleID}
fi
tmpDir=$outDir/tmp.$SLURM_JOB_ID
if [ ! -d "$tmpDir" ]; then
    mkdir -p $tmpDir
fi

# Load modules
for mod in "${modList[@]}"; do
    module load $mod
done

# Revert BAMs to fastq
echo "## INFO: Processing sample: ${sampleID}" # helps with troubleshooting array jobs
samtools sort -l 0 -m 4G -n -@${nCores} -T$tmpDir --reference ${genomeBuild} ${bamFile[SLURM_ARRAY_TASK_ID]} |\
samtools fastq -1 $outDir/${sampleID}/${sampleID}.reads_R1.fastq.gz -2 $outDir/${sampleID}/${sampleID}.reads_R2.fastq.gz -0 /dev/null -s $outDir/${sampleID}/${sampleID}.reads_U1.fastq.gz -n -@${nCores} -

# Clean up
rm -r $tmpDir
