# Metabolomics-Data-Integration

## Background and Approach
Untargeted metabolomics data is challenging to interpret since it is typically comprised of thousands of known and unknown metabolites. An additional challenge is in integrating this high-dimensional data with other high-dimensional data (e.g., microbiome, transcriptome). Thus, I applied a data reduction technique, weighted co-network analysis, to group metabolites that are highly correlated with one another into modules. I applied the [WGCNA package in R](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) to metabolomics data (after log-transformation and filtering). After grouping metabolites into modules, an Eigennode is calculated for each module for each biological replicate, which allows for modeling as a continuous variable in statistical models. 

## Integration with Phenotyping Data
The R code provided here demonstrates the use of modeling metabolite module Eigenvalues in linear mixed effects regression (LMER) models to investigate the relationships between the gut microbiome and phenotypes and genotypes longitudinally. An additional R script is provided for generating heatmaps depicting these associations. Below is an example:

[heatmap]()


## Integration with Microbiome Data
An approach for the integration of metabolomics data with microbiome data is in the [Microbiome-Metabolome-Integration Repository.](https://github.com/karineier/Microbiome-Metabolome-Integration) 

## Application to Biomedical Research
I applied these analyses to investigate molecular signatures of disease progression in a mouse model of Rett syndrome, a neurodevelopmental disorder that features metabolic and gastrointestinal dysfunction. The findings are published in [Communications Biology.](https://rdcu.be/cDkCI)

