#!/bin/bash

# Example usage:
# PREFIX=ExomesMarch2014 INDIR=/data/neurogenetics/gVcfDumpingGround/Exomes OUTDIR=/data/neurogenetics/variants/gvcf sbatch /data/neurogenetics/git/PhoenixScripts/shared/scripts/GATKGenoRecal.split.q -J GATKGeno

#SBATCH -J GATKGeno

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=3-00:00
#SBATCH --mem=80GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# load modules
module load Java/1.8.0_121
module load HTSlib/1.3.1-GCC-5.3.0-binutils-2.25
module load R/3.4.2-foss-2016b 

# run the executable
/data/neurogenetics/git/PhoenixScripts/shared/scripts/GATKv3.x.VCFGenotypeRecal.split.HPC.pipeline.sh -p $PREFIX -i $INDIR -o $OUTDIR
