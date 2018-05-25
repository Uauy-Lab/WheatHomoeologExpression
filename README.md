# WheatHomoeologExpression

# Data 

All the input data for the scripts described can be found [here](https://opendata.earlham.ac.uk/wheat/under_license/toronto/).


## Methods for Transcriptome paper. 

This repository contains the code as used for the following sections

 <span style="display: none;">
### Mapping and reference

</span>

### Homoeolog specific mapping by kallisto in nulli-tetrasomic lines

* [plot heatmaps for nullitetra, all genes (not just triads)](scripts/nullitetra/1_plot_nullitetra_heatmap.R)
* [get expression for nullitetras 1:1:1 triads](scripts/nullitetra/2_get_expression_nullitetra_1_1_1_triads.R)
* [add dominance category to nullitetra 1:1:1 triad expression](scripts/nullitetra/3_add_dominance_category_to_1_1_1_triads.R)
* [plot and stats for A vs B vs D for nullitetras and dominance categories](scripts/nullitetra/4_boxplot_1_1_1_triads_A_vs_D_vs_D_and_dominance_final_1tpm_threshold.R)

 <span style="display: none;">
### Theoretical read mapping accuracy between homoeologs using SNP distributions

</span>

### Genome of origin effect on gene expression

* [plot clustering of intermediate tissue expression, colour by genome of origin](scripts/genome_of_origin_effect/hclust_ABD_UPLOAD.r)

### Expression complexity
Details in Jupyter [notebook](05.%20Cumulative%20expression.ipynb)  

### Expressed genes within the Azhurnaya developmental time course (209 samples) and Chinese Spring no stress (123 samples).
    
Details in Jupyter [notebook](01.%20Prepare%20TPM%20table.ipynb)  


### Differential expression (Azhurnaya time course)

* [differential expression all tissues compared to anther](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_anther.R)
* [differential expression all tissues compared to awns](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_awns.R)
* [differential expression all tissues compared to embryo](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_embryo.R)
* [differential expression all tissues compared to endosperm](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_endosperm.R)
* [differential expression all tissues compared to flag leaf blade](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_flag_leaf_blade.R)
* [differential expression all tissues compared to flag leaf sheath](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_flag_leaf_sheath.R)
* [differential expression all tissues compared to glumes](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_glumes.R)
* [differential expression all tissues compared to grain hard](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_grain_hard.R)
* [differential expression all tissues compared to grain milk](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_grain_milk.R)
* [differential expression all tissues compared to internode](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_internode.R)
* [differential expression all tissues compared to leaf blades excl flag](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_leaf_blades_excl_flag.R)
* [differential expression all tissues compared to leaf ligule](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_leaf_ligule.R)
* [differential expression all tissues compared to leaf sheaths excl flag](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_leaf_sheaths_excl_flag.R)
* [differential expression all tissues compared to peduncle](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_peduncle.R)
* [differential expression all tissues compared to root apical meristem](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_root_apical_meristem.R)
* [differential expression all tissues compared to seedling aerial tissues](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_seedling_aerial_tissues.R)
* [differential expression all tissues compared to shoot apical meristem](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_shoot_apical_meristem.R)
* [differential expression all tissues compared to shoot axis](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_shoot_axis.R)
* [differential expression all tissues compared to spike](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_spike.R)
* [differential expression all tissues compared to spikelets](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_spikelets.R)
* [differential expression all tissues compared to stigma ovary](scripts/Tissue_differential_expression/1_DE_intermediate_tissues_script_for_cluster_one_tissue_at_a_time_stigma_ovary.R)
* [after differential expression analysis filter differentially expressed genes for padj < 0.05](scripts/Tissue_differential_expression/2_filter_DE_padj0.05.R)
* [after filtering padj < 0.05 summarise differentially expressed genes ](scripts/Tissue_differential_expression/3_summarise_DE_intermediate_tissues.R)


### eFP browser 

The browser is available [here](http://bar.utoronto.ca/efp_wheat/cgi-bin/efpWeb.cgi)

### Relative expression levels of the A, B, and D genome homoeologs across triads

Details in Jupyter [notebook](02.%20Calculate%20triad%20category.ipynb)  

### Definition of homoeolog expression bias categories

Details in Jupyter [notebook](02.%20Calculate%20triad%20category.ipynb#Definition-of-homoeolog-expression-bias-categories)  

### Analysis of effect on polyploidy on homoeolog expression bias
Details in Jupyter [notebook](02.%20Classify%20genes%20by%20movement) 

### Variation in homoeolog expression bias across tissues (Stable/Dynamic triads).
Details in Jupyter [notebook](04.%20Synthetic%20wheat%20calculate%20TPMs.ipynb)

### Gene Ontology and Plant Ontology term enrichment
Ricardo

### CDS, protein and promoter analysis for triads
The scripts used for this section are in various repositories: 

  * [Extract promoter regions](https://github.com/lucventurini/mikado/blob/f47aa63/util/extract_promoter_regions.py)
  * [Pairwise alignment of triads](https://github.com/TGAC/bioruby-polyploid-tools/blob/master/bin/blast_triads.rb)

### Transposable Element (TE) modeling using CLARITE

Clarite is available [here](https://github.com/jdaron/CLARI-TE)

### TE presence in gene promoters

* [Calculate sliding window of Transposable Element density upstream of genes](scripts/TEs/calculate_sliding_windows_TEs.r)
* [Calculate average lengths and distances from start site of Transposable Elements including statistical tests](scripts/TEs/get_average_length_of_TEs_and_distance_to_ATG.r)
* [Script for analysis of TE presence/absence based on homoeolog expression bias categories.](scripts/TEs/retro_presenceabsence_subsetting_analysis_dominancecategories_with1500_UPLOAD.R)
* [Script for analysis of TE presence/absence based on homoeolog expression bias variation categories.](scripts/TEs/retro_presenceabsence_subsetting_analysis_movementcategories_with1500_UPLOAD.R)


### Enrichment of TE families in gene promoters

* [Script for analysis of TE family enrichment based on homoeolog expression bias variation categories.](scripts/TEs/TE_Family_Comparison_10-80-10_Movement_UPLOAD.R)

### TE density in gene promoters

* [Plot sliding window of Transposable Element density including statistical tests](scripts/TEs/plot_final_TE_density_graphs.r)

### Transcription factor binding site (TFBS) identification

* [Analyse TFBS categories according to triad movement including statistical tests](scripts/TFBS/analyse_TF_motif_categories_1500bp.r)
* [Filter unique FIMO hits only - step 1](scripts/TFBS/create_FIMO_unique_hits_tables_1500_slurm.sh)
* [Generate TFBS categories based on occurrence in triads](scripts/TFBS/generate_TF_motif_groupings_1500bp.r)
* [Filter unique FIMO hits only - step 2](scripts/TFBS/generate_unique_TF_motifs_hits_1500bp.r)
* [Run FIMO to identify TFBS in promoters](scripts/TFBS/IWGSCv1.0_UTR.HC_1500bp_upstream_FIMO.sh)


### Ka/Ks analysis

 * [MAFFT alignment](https://github.com/TGAC/bioruby-polyploid- tools/blob/master/bin/mafft_triads.rb.)
 * [Script to obtain Chinese Spring No-Stress analysis of Ka/Ks ratios for syntenic 1:1:1 triads.](scripts/KaKs/kaks_chinesespring_UPLOAD.R)
 * [Script to obtain Azhurnaya development analysis of Ka/Ks ratios for syntenic 1:1:1 triads.](scripts/KaKs/kaks_developmentset_UPLOAD.R)
 * [Script to obtain Chinese Spring No-Stress analysis of Ka/Ks ratios for non-syntenic 1:1:1 triads.](scripts/KaKs/kaks_nonsyntenictriads_chinesespring_UPLOAD.R)
 * [Script to obtain Azhurnaya development analysis of Ka/Ks ratios for non-syntenic 1:1:1 triads.](scripts/KaKs/kaks_nonsyntenictriads_development_UPLOAD.R)


### WGCNA Network construction 

* [WGCNA combine studies filter normalise  samples for abiotic.R](scripts/WGCNA_network_construction/1_WGCNA_combine_studies_filter_normalise_cluster_abiotic.R)
* [WGCNA combine studies filter normalise  samples for disease.R](scripts/WGCNA_network_construction/1_WGCNA_combine_studies_filter_normalise_cluster_disease.R)
* [WGCNA combine studies filter normalise  samples for grain.R](scripts/WGCNA_network_construction/1_WGCNA_combine_studies_filter_normalise_cluster_grain.R)
* [WGCNA combine studies filter normalise  samples for leaf.R](scripts/WGCNA_network_construction/1_WGCNA_combine_studies_filter_normalise_cluster_leaf.R)
* [WGCNA combine studies filter normalise  samples for root.R](scripts/WGCNA_network_construction/1_WGCNA_combine_studies_filter_normalise_cluster_root.R)
* [WGCNA combine studies filter normalise  samples for spike.R](scripts/WGCNA_network_construction/1_WGCNA_combine_studies_filter_normalise_cluster_spike.R)
* [WGCNA determine threshold for abiotic.R](scripts/WGCNA_network_construction/2_WGCNA_cluster_thresholding_abiotic.R)
* [WGCNA determine threshold for disease.R](scripts/WGCNA_network_construction/2_WGCNA_cluster_thresholding_disease.R)
* [WGCNA determine threshold for grain.R](scripts/WGCNA_network_construction/2_WGCNA_cluster_thresholding_grain.R)
* [WGCNA determine threshold for leaf.R](scripts/WGCNA_network_construction/2_WGCNA_cluster_thresholding_leaf.R)
* [WGCNA determine threshold for root.R](scripts/WGCNA_network_construction/2_WGCNA_cluster_thresholding_root.R)
* [WGCNA determine threshold for spike.R](scripts/WGCNA_network_construction/2_WGCNA_cluster_thresholding_spike.R)
* [WGCNA run WGCNA abiotic.R](scripts/WGCNA_network_construction/3_WGCNA_cluster_co-expression_abiotic.R)
* [WGCNA run WGCNA disease.R](scripts/WGCNA_network_construction/3_WGCNA_cluster_co-expression_disease.R)
* [WGCNA run WGCNA grain.R](scripts/WGCNA_network_construction/3_WGCNA_cluster_co-expression_grain.R)
* [WGCNA run WGCNA leaf.R](scripts/WGCNA_network_construction/3_WGCNA_cluster_co-expression_leaf.R)
* [WGCNA run WGCNA root.R](scripts/WGCNA_network_construction/3_WGCNA_cluster_co-expression_root.R)
* [WGCNA run WGCNA spike.R](scripts/WGCNA_network_construction/3_WGCNA_cluster_co-expression_spike.R)


### Defining same, similar and divergent expression patterns of triads

* [calculate the distance between homoeologs](scripts/same_similar_divergent_homoeolog_expression/1_calculate_distance_between_homoeologs_all_tissues.R)
* [plot graphs of distances between homoeologs](scripts/same_similar_divergent_homoeolog_expression/2_plot_graphs_distance_all_tissues_triads_for_manuscript.R)

### Module overlaps

* [compares whether modules in different tissue networks have significant overlap in terms of genes in each module](scripts/compare_tissue_modules/compare_modules.R)


### Gene tree for MADS_II in root module 61
* [blast MADSII proteins against other wheat proteins](scripts/MADSII_root/blast_MADSII_root_module61.sh)


### Correlation to stress status

* [calculate stress correlations for abiotic](scripts/stress_correlations/module_trait_relationships_padj_abiotic.R)
* [calculate stress correlations for disease](scripts/stress_correlations/module_trait_relationships_padj_disease_merge_same_stress.R)


### Genie3

* [Results from Genie3](https://doi.ipk-gatersleben.de/DOI/53148abd-26a1-4ede-802b-c2635af6a725/0dd8224a-34fc-4e3b-8ab8-883d07e52bd2/2/1847940088)

### Identifying high connected hub genes

* [find hub genes in abiotic network](scripts/hub_genes/find_hub_genes_abiotic.R)
* [find hub genes in disease network](scripts/hub_genes/find_hub_genes_disease.R)




