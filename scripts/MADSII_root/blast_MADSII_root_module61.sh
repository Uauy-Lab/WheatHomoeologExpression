#!/bin/bash
#
# SLURM batch script to launch BLAST
#
#SBATCH -p nbi-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --mem 3000 # memory pool for all cores
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH -o /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/tissues/transporters/slurm_output/extract_CSS_genes.%N.%j.out # STDOUT
#SBATCH -e /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/tissues/transporters/slurm_output/extract_CSS_genes.%N.%j.err # STDERR
#SBATCH -J blast
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=philippa.borrill\@jic.ac.uk # send-to address

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/WGAv1.0/annotation

source blast+-2.2.28
makeblastdb -in iwgsc_refseqv1.0_HighConf_PROTEIN_2017Mar13.fa -dbtype prot -parse_seqids

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/expression_browser/WGA/WGCNA/tissues/root/maxP0.05/gene_tree_module61_MADS_II

blastp -db /nbi/Research-Groups/NBI/Cristobal-Uauy/WGAv1.0/annotation/iwgsc_refseqv1.0_HighConf_PROTEIN_2017Mar13.fa -query wheat_protein_MADSII.txt -num_threads 1  -max_target_seqs 30 -outfmt 6 -out wheat_MADSII_module61_vs_WGAv1.0_peptide_max30targets.txt

# run this after making a list of the top 30 hits from wheat
#blastdbcmd -db /nbi/Research-Groups/NBI/Cristobal-Uauy/WGAv1.0/annotation/iwgsc_refseqv1.0_HighConf_PROTEIN_2017Mar13.fa -entry_batch list_of_wheat_MADSII_top_30_hits.txt -out wheat_MADSII_protein_top30.fa
