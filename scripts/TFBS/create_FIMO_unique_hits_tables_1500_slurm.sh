#!/bin/bash -e
#SBATCH -p RG-Cristobal-Uauy
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem 5000
#SBATCH -o /nbi/group-data/ifs/NBI/Cristobal-Uauy/Jemima/companion_paper/slurm_out/FIMO_unique_hits.%N.%j.out
#SBATCH -e /nbi/group-data/ifs/NBI/Cristobal-Uauy/Jemima/companion_paper/slurm_out/FIMO_unique_hits.%N.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jemima.brinton@jic.ac.uk

#generate file with unique binding sites only
FIMO_RESULTS_PATH='/nbi/group-data/ifs/NBI/Cristobal-Uauy/Jemima/companion_paper/out_2018_03_12_12_03/fimo_results/'

cd $FIMO_RESULTS_PATH

FIMO_RESULTS=($(ls $FIMO_RESULTS_PATH))

for RESULTS_DIR in ${FIMO_RESULTS[*]}
	do
	FIMO_RESULTS=$RESULTS_DIR'/fimo.txt'
	
	#create an analysis directory
	OUT_DIR=$RESULTS_DIR'/unique_hits_analysis/'
	
	mkdir -p $OUT_DIR
	
	#generate a file that just contains the unique hits
	
	cut -f 1,2,3 $FIMO_RESULTS | uniq -u > $OUT_DIR'fimo_unique.tsv'
	
	#now to use generate_unique_TF_motifs_hits.r to furhter simplify the output
	done

