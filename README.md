# Supplementary repository to replicate the coLCNEC analysis in the paper: Bronco-Pulmonary Combined large-cell neuroendocrine carcinomas (Co-LCNEC): more than a hybrid neoplasm.

********************************
***CONTENT OF THIS REPOSITORY***
********************************
- R scripts to replicate the methodological section included in the paper about gene expression data
- Description of each file in the repository

********************************
***DESCRIPTION OF WORKFLOW***
********************************

![Test Image 8](https://github.com/LucaGiudice/Supplementary-coLCNEC/blob/main/data/workflow.png)

********************************
***DESCRIPTION OF EACH FILE AND SUBDIRECTORY***
********************************

- ***/op3_dga.R*** Find core gene sets for each starting sample group and find the clustering setting which minimaze the number of intruders in each final cluster

- ***/op4_CL.R*** Perform the final clustering

- ***/op4_data_prep.R*** Prepare the data for the pairwise comparisons between the clusters
   
- ***/op4x_artefatti.R*** Find the outlier genes to exclude from the results of the analysis

- ***/op5_DE_analysis.R*** Perfrom the DE analysis of all the defined contrasts

- ***/op6_GSEA.R*** Perform gene set enrichment analysis with GSEA, SPIA, ORA and Camera

- ***/op7_report.R*** Produce the tables and files of the results

- ***/op8_table_genes.R*** Produce the summerized table of the genes for the paper

- ***/op8_table_gsea.R*** Produce the summerized table of the pathways for the paper

- ***/op9_gene_overexp.R*** Perform an analysis to get overexpressed genes between clusters (excl1 only for one cluster, excl2 for a pair of clusters)

- ***/op9_pathways_overexp.R*** Perform an analysis to get overexpressed pathways using the ssGSEA measure (c6 only for c6 pathways, h with only the Hallmark pathways)

- ***/op10_final_gene_heatmap.R*** Produce the final heatmap of the samples in study

- ***/op11_final_gene_list.R*** Produce the final gene list

- ***/op11_final_pathways_list.R*** Produce the final pathway list

- ***/op12_umap.R*** Produce the UMAP of the samples, classes and clusters
