#!/bin/bash

# Example usage:
# <variables> sbatch BWA-GATKHC.forGenomes.q

#SBATCH -J GenomeBWA-GATKHC
#SBATCH -o /hpcfs/users/%u/launch/slurm-%j.out

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=1:00:00 # change this to 3 days for real set
#SBATCH --mem=24GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=alison.gardner@adelaide.edu.au

# load modules
module load arch/haswell
module load BWA/0.7.17-foss-2016b
module load Java/1.8.0_121
module load HTSlib/1.10.2-foss-2016b
module load SAMtools/1.10-foss-2016b
### module load GATK 3.7
### module load picard/2.6.0 or higher

# run the executable
/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/GATK4/BWA-GATKHC.Nov2020_hs37_Phoenix.sh -p 12200 -s /hpcfs/groups/phoenix-hpc-neurogenetics/sequences/Illumina/exome/CP-Kruer/12200/concat/ -o /hpcfs/groups/phoenix-hpc-neurogenetics/sequences/Illumina/exome/CP-Kruer/12200/concat/ -c /hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/GATK4/BWA-GATKHC.GRCh37_full_analysis_set_phoenix.cfg -S 12200 -L IlluminaGenome -I 12200
