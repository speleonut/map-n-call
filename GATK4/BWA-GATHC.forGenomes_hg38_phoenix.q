#!/bin/bash

# Example usage:
# <variables> sbatch BWA-GATKHC.forGenomes.q

#SBATCH -J GenomeBWA-GATKHC
#SBATCH -o /fast/users/a1149120/launch/slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=2:00:00 # change this to 3 days for real set
#SBATCH --mem=24GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=alison.gardner@adelaide.edu.au

# load modules
module load BWA/0.7.15-foss-2017a
module load Java/1.8.0_121
module load HTSlib/1.9-foss-2016b
module load SAMtools/1.9-foss-2016b
### module load GATK 3.7
### module load picard/2.6.0 or higher

# run the executable
/data/neurogenetics/git/PhoenixScripts/shared/scripts/BWA-GATKHC.Apr2020_hg38_Phoenix.sh -p TESTX -s /fast/users/a1149120/Test_fastqs/ -o /fast/users/a1149120/Test_fastqs -c /data/neurogenetics/git/PhoenixScripts/shared/scripts/BWA-GATKHC.GRCh38_full_analysis_set_phoenix.cfg -S TestX -L IlluminaGenome -I TestX
