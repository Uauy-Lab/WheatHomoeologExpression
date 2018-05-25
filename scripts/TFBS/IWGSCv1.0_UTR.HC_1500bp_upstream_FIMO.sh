#!/bin/bash -e
#SBATCH -p RG-Cristobal-Uauy
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem 100000
#SBATCH -o ./out_2018_03_12_12_03/run_logs/FIMO.%N.%j.out
#SBATCH -e ./out_2018_03_12_12_03/run_logs/FIMO.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jemima.brinton@jic.ac.uk

cd /nbi/Research-Groups/NBI/Cristobal-Uauy/Jemima/companion_paper/

source meme-4.11.4

#create background file

fasta-get-markov /nbi/Research-Groups/NBI/Cristobal-Uauy/WGAv1.0/annotation/upstream_cds_fasta_files/fullset_consensus_firstCDS_USETHESE/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa out_2018_03_12_12_03/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa_background

#run fimo_results
fimo --bgfile out_2018_03_12_12_03/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa_background --thresh 1e-4 --motif-pseudo 0.00000001 --max-stored-scores 1000000 --o out_2018_03_12_12_03/fimo_results/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa /nbi/Research-Groups/NBI/Cristobal-Uauy/Jemima/5A_transcriptomics_manuscript/outside_region_promoter_analysis/FIMO/plantPAN2.0/PlantPAN2_TFBS_wheat.meme /nbi/Research-Groups/NBI/Cristobal-Uauy/WGAv1.0/annotation/upstream_cds_fasta_files/fullset_consensus_firstCDS_USETHESE/IWGSCv1.0_UTR.HC.firstCDS_consensus._1500bp_upstream_triads_only.fa
