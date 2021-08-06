
#### using linear mixed effects regression to model longitudinal associations between WGCNA modules and traits ####

setwd("/Users/karineier/Documents/Mecp2/WGCNA_Metabolomics/LMER_models")

### Loading libraries needed ###

library("WGCNA")
options(stringsAsFactors=FALSE)
library("ggplot2")
library("lme4")
library("lmerTest")
library("tidyverse")
library("openxlsx")

longitudinalModels = function(sex) {
  
  sexname.0 = ifelse(sex=="F", "Females", "Males")
  sexname.1 = ifelse(sex=="F", "Female", "Male")
  
  # Load WGCNA data from part 0 and part 1 #
  
  load(glue::glue("/Users/karineier/Documents/Mecp2/WGCNA_Metabolomics/{sexname.0}/{sexname.1}FecalMetabolomics-01-dataInput.RData"))
  load(glue::glue("/Users/karineier/Documents/Mecp2/WGCNA_Metabolomics/{sexname.0}/{sexname.0}-02-networkConstruction-auto.RData"))
  
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0) 
  
  MEs.1 = MEs[
    with(MEs, order(rownames(MEs))),
  ]
  
  traits = read.csv("/Users/karineier/Documents/Mecp2/WGCNA_Metabolomics/Metabolomics_Traits.csv")
  traits = subset(traits, Sex_cat==sex)
  
  traits = traits[
    with(traits, order(traits$Sample_ID)),
  ]
  
  traits$Neuro_phenotype = traits$Neuro_phenotype +1 ### Adding 1 to Neurophenotyping scores to avoid issues with scores of 0
  
  alldata = cbind(traits, MEs.1)
  alldata$Mouse_ID = as.factor(alldata$Mouse_ID)
  
  head(alldata)
  
  traitcor = cor(MEs.1, traits[,-c(1:4, 6, 7)], use="p", method="spearman") 
  
  openxlsx::write.xlsx(traitcor, file=glue::glue("traitcors_{sex}.xlsx"))
  
  ### setting up models ###
  
  allModels = function(traitname) {
    
    fitlist = as.list(1:no.MEs)
    names(fitlist) <- MEnames
    
    p.values = data.frame(matrix(ncol=1, nrow=no.MEs))
    rownames(p.values) = MEnames
    colnames(p.values) = "p.value"
    
    betas = data.frame(matrix(ncol=1, nrow=no.MEs))
    rownames(betas) = MEnames
    colnames(betas) = "beta"
    
    
    for(i in MEnames){
      
      # print status
      print(paste("Running entity:", i, "which is", which(MEnames==i), "out of", no.MEs))
      
      #create temporary data matrix and model formula
      
      fml <- as.formula( paste(i, "~", traitname, "+ (1|Mouse_ID)"))
      
      #assign fit to list by name
      fitlist[[i]] <- lmer(fml, data=alldata)
      p.values[i,1] <- summary(fitlist[[i]])$coefficients[2,5]
      betas[i,1] <- summary(fitlist[[i]])$coefficients[2,1]
      
      
    }
    
    results=cbind(p.values, betas)
    
    return(results)
    
  }
  
  list = lapply(traitnames, allModels)
  names(list) = traitnames
  
  all = do.call(cbind.data.frame, list)
  
  number = ncol(all)
  
  odd <- seq(1, number, 2)
  even <- seq(2, number, 2)
  
  all.pvals = all[,odd]
  all.betas = all[,even]
  
  all.pvals.mat = as.matrix(all.pvals)
  all.FDRs = matrix(p.adjust(as.vector(all.pvals.mat), method="fdr"), ncol=ncol(all.pvals)) 
  all.FDRs = as.data.frame(all.FDRs)
  colnames(all.FDRs) = traitnames
  all.FDRs$ME = MEnames
  
  colnames(all.pvals) = traitnames
  all.pvals$ME = MEnames
  
  colnames(all.betas) = traitnames
  all.betas$ME = MEnames
  
  openxlsx::write.xlsx(all.pvals, file=glue::glue("P.values_{sex}.xlsx"))
  openxlsx::write.xlsx(all.FDRs, file=glue::glue("FDRs_{sex}.xlsx"))
  openxlsx::write.xlsx(all.betas, file=glue::glue("Betas_{sex}.xlsx"))
  
  save(all.pvals, all.FDRs, all.betas, traitcor, file=glue::glue("LMER_models_traits_v_WGCNA_{sex}.RData"))
  
}

lapply(c("F", "M"), longitudinalModels)




