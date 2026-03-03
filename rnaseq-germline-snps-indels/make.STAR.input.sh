#!/bin/bash
#SBATCH -J makeSTAR
#SBATCH -o /hpcfs/users/%u/log/makeSTARinput.slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1                            # number of nodes (due to the nature of sequential processing, here uses single node)
#SBATCH -n 1                            # number of cores (here uses 12)
#SBATCH --time=00:10:00                 # time allocation, which has the format (D-HH:MM)
#SBATCH --mem=1G                        # memory pool for all cores (here set to 32 GB)

# Notification configuration
#SBATCH --mail-type=END                 # Type of email notifications will be sent (here set to END, which means an email will be sent when the job is done)
#SBATCH --mail-type=FAIL                # Type of email notifications will be sent (here set to FAIL, which means an email will be sent when the job is fail to complete)
#SBATCH --mail-user=%u@adelaide.edu.au  # Email to which notification will be sent

# A short script to make up input for STAR.
# This script is optimised to work with the GATK4 Universal RNAseq Germline SNPs and Indels pipeline but can be used for other purposes if needed. YMMV
# The script expects each sample to have a folder with R1 and R2 fastq.gz files.
# Usage: sbatch make.STAR.input.sh /path/to/sequence/file/folders
# 
# Forked from https://github.com/speleonut/RNASeq/blob/master/hisat2.Illumina/makeHISAT2input.sh
# Mark Corbett 03/03/2026

# Check to see if script has the correct input
if [ -z "$1" ]; then
    echo "## ERROR: You didn't specify the path to your sequence files.
    Usage: $0 /path/to/sequence/file/folders"
    exit 1
fi

if [ ! -d "$seqPath" ]; then
    echo "Usage: $0 /path/to/sequence/file/folders
    ## ERROR: $seqPath is not a directory"
    exit 1
fi

# Fetch sample ID from the folder names.

ls -d */ | cut -f1 -d"/" | sort | uniq > $seqPath/tmp.make.STAR.input.SampleID.txt # Temporary filename

while read i; do
	find ${seqPath}/${i}/ -name "*.gz" > $i.files.txt
	grep R1.f $i.files.txt > R1.tmp.txt
	grep R2.f $i.files.txt > R2.tmp.txt
	rm $i.files.txt
	# for this next bit the code comes from
	# http://unix.stackexchange.com/questions/114943/can-sed-replace-new-line-characters
	# it works!
	sed ':a;N;$!ba;s/\n/,/g' R1.tmp.txt >> R1.tmp.list.txt
	sed ':a;N;$!ba;s/\n/,/g' R2.tmp.txt >> R2.tmp.list.txt
	rm R1.tmp.txt R2.tmp.txt
done < $seqPath/tmp.make.STAR.input.SampleID.txt

paste $seqPath/tmp.make.STAR.input.SampleID.txt R1.tmp.list.txt R2.tmp.list.txt > $seqPath/STAR.input.list.txt
rm R1.tmp.list.txt R2.tmp.list.txt $seqPath/tmp.make.STAR.input.SampleID.txt
