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

We prepared the data composed of 20,815 genes following the workflow defined by Law et al.[30], performing a quality control and removing the sample’s differences based on how meta information drove the count values. Genes not expressed at a biologically meaningful level in any condition were filtered out to increase the reliability of the mean-variance relationship. In this dataset, the median library size was about 5.8 million, so the filterByExpr function kept genes that had a log-CPM of 1.55 or more in at least three samples. The differences between samples due to the depth of sequencing were then removed, the data were normalised using the trimmed mean of M-values (TMM)[31] method and the count values transformed to log-CPM. We proceeded to check for outsider samples. We assessed how much the single sample medians were distant from the overall one using the z-score. We filtered out 2 samples having a score greater than 3 indicating a distribution of expression values too different from the remaining cases. The final dataset was composed of 128 samples per 18616 genes and was manipulated in two different ways during the downstream analysis. In one case, it was batch corrected, used to find and to pathway-level characterize sample clusters. In the second case, it was not corrected but used directly to get cluster specific enriched genes with the meta information integrated as covariates.
For the first scenario, we treated the dataset following the pipeline included in the package proBatch[32] to check if location, material and batch were biasing the expression values and influencing the grouping of the samples. We performed agglomerative hierarchical clustering Agnes[33, 34] to diagnose for batch effects and to evaluate to what extent technical variances still existed in the normalized data. We estimated the method (average, single, complete, ward) which maximized the clustering coefficient, we determined the optimal number of clusters using the average silhouette width and we applied the clustering algorithm. The resulting clusters have been compared using the adjusted Rand index[35] to the samples separated respectively by the material, location and batch information. The batch variable got an index of 0.25 higher than the other information and together with the visualization of the samples by the projection of the gene values on two principal components highlighted that the batch information was driving the grouping of the samples. We then corrected the batch affect using ComBat, described in Johnson et al.[36]. The normalized and corrected dataset is used for a semi-supervised consensus clustering to unveil the similarities inside and between the histological classes. For this operation, we kept only the genes which had not overlapping expression distributions between the histological groups to some extent for modelling a degree of tolerance. The gene distribution composed of its expression values in a group has been summarized with an interval defined by two limits equal to the 0.25 (lower bound) and the 0.75 (upper bound) quantiles. A gene has been selected and defined core for a group when it had a distribution not overlapping with its intervals in the other groups (at least the 60% of the groups). This allowed us to catch the differences and similarities between samples represented by their own almost-core genes without leading the clustering algorithm to provide overfitted results (i.e. clusters perfectly separated due to the starting histological classes). The consensus clustering has been performed with Cola[37]. The method tests six partitioning methods combined to four different feature selection strategies for each possible k number of clusters. After the tests, it provides four scores (Silhouette score, PAC score, Concordance and Jaccard index) that allow to select the best combination. In our case, the method Pam with the standard deviation as approach to retrieve the most informative genes and k equal to ten got the best clustering after that has been tested in 50 runs with resample each time of the 80% of the genes. After having unveiled the real groups of samples, we removed outlier genes for not including them in the downstream analysis. We estimated the coefficients of skewness, kurtosis and linearity together with the position and value of the highest peak in the distributions assumed by each gene in the clusters. We removed the genes that in at least three clusters had a fat-tailed distribution with a kurtosis coefficient greater than 3, a non-linear distribution with a coefficient lower than 1, a peak in the tail above the 0.95 quantile and a peak twice higher the median. These parameters have been found by performing a tuning operation which criteria was to maximize the removal of few manually selected outliers characterized by a low variation in the expression values but a peak in a number of samples lower than the 5% in multiple clusters. Next, we moved to the gene and pathway analysis. 
We created a design matrix with the location and material as covariates for each pair of clusters. We computed differential expressed genes following limma workflow. We intersected the results of all the pairwise comparisons and determined how much a gene was frequently appearing as DE for a specific cluster against every other and poorly appearing in the comparisons not including the group. This step allowed us to get DE cluster-specific genes. We finished the analysis at the gene-level determining the overexpressed genes which were exclusive of one cluster. In other words, the genes which in a group had the distribution described by the 0.3 quantile greater than its 0.8 quantile in the other clusters. We replicated the same step considering the exclusive overexpression of a gene also in a pair of clusters against all the others. In this last case, the distribution of the gene could overlap completely between the clusters of the pair.
Next, we moved to analyse the samples at the pathway-level. We downloaded c6 and H pathways from MSigDB[38, 39] and determined the cluster-specific enriched gene sets using the normalized and batch corrected count matrix. We applied GSEA using GAGE R package[40] between clusters to get pairwise significant up and down regulated pathways. While, we used an approach based on the ssGSEA score[38] for determining the biological processes differently enriched between all the clusters. The ssGSEA score determines how much the genes in a particular set are co-ordinately up- or down-regulated within a specific sample. We assessed the ssGSEA score for each pair of sample and gene set. We represented each pathway in a cluster with the mean score got by its members. We performed a z-score normalization of the pathway scores in the clusters. We ranked the biological processes for each cluster based on the Euclidian distance between its enrichment score and the mean of the other groups. We then selected the top pathways in rank list to characterize the clusters. We downloaded the KEGG signaling pathways[41, 42] and performed SPIA[43] to assess which one was perturbed by the set of DE genes which were significant in each designed contrast. All the analysis and results can be reproduced with the R scripts shared at the following github page: https://github.com/LucaGiudice/Supplementary-coLCNEC

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
