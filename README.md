# map-n-call
Scripts for mapping DNA sequencing and calling variants

## Mapping Illumina short reads to the genome and calling variants
The following is an adaptation of the BWA - GATK4 "best practices" guidelines (or at least somewhat like it) for the University of Adelaide Phoenix / HPC1 facility.

## Getting started
It's best to use our Google Doc to get started:

https://docs.google.com/document/d/1TbvCCqPcnCx9Xl9MZPZgQquyuG6FG43shAobZbsFB8I/edit#heading=h.l554c6o77bhn

There are two pipelines, each composed of dependent SLURM jobs. The first covers mapping with BWA to generating a gVCF with GATK4 haplotype caller for a single sample. This pipeline is launched using the `GATK4/BWA-GATKHC.Universal_Launcher.Phoenix.sh` script.

The second performs genotyping over multiple samples and classifies variants using VQSR. This pipeline is launched using the `GATK4/GATK.GenoRecal_Universal_Launcher_phoenix.sh` script

## Pipeline defaults
If you don't tell them otherwise these scripts will automatically use the following defaults:

Genome build:  hs38DH

BAM file output location: `/hpcfs/users/${USER}/BWA-GATKHC/${Sample}`

gVCF file output location: `/hpcfs/groups/phoenix-hpc-neurogenetics/variants/gVCF/${BUILD}`

VCF file output location: `/hpcfs/groups/phoenix-hpc-neurogenetics/variants/vcf/${BUILD}`

SLURM log file location: `/hpcfs/users/${USER}/log` (note the change from "launch" to "log" compared with previous versions of this script)

GATK4 log file locations: These will end up with your BAMs in the first pipeline and the final VCF in the second pipeline

## Making new config files
All config files have been moved to the `configs/` directory.  Only files with the file suffix _phoenix.config work with the GATK4 pipeline.  To make a new config please use the supplied `configs/BWA-GATKHC.TEMPLATE_phoenix.cfg` file to guide you *and* update (or request an update) the `test_genome_build` function in the `GATK4/BWA-GATKHC.Universal_Launcher.Phoenix.sh` script if required.

Happy mapping & calling!