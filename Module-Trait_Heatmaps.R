
####### Module-Trait Association Heatmaps Using FDRs from LMER models ###########

setwd("/Users/karineier/Documents/Mecp2/WGCNA_Metabolomics/LMER_models")

library("WGCNA")
library("tidyr")

traitHeatmaps = function(sex) {
  load(glue::glue("LMER_models_traits_v_WGCNA_{sex}.RData"))
  
  rownames(all.FDRs) = all.FDRs$ME
  all.FDRs = as.matrix(all.FDRs[,-c(ncol(all.FDRs))])
  
  textMatrix.traits = paste(ifelse((signif(all.FDRs, 1))<0.05, "*", (ifelse((signif(all.FDRs, 1))<0.1, "^", ""))), sep="")
  textMatrix.traits = matrix(textMatrix.traits, ncol=ncol(all.FDRs))
  rownames(traitcor) = gsub("ME", "", rownames(traitcor))
  
  tiff(glue::glue("Trait_Heatmap_Full_{sex}.tiff"), res=400, height=5, width=6.5, units="in")
  map1 = labeledHeatmap(Matrix = traitcor,
                        xLabels = colnames(traitcor),
                        yLabels = rownames(traitcor),
                        ySymbols = rownames(traitcor),
                        colorLabels = FALSE,
                        colors = blueWhiteRed(50),
                        textMatrix = textMatrix.traits,
                        setStdMargins = FALSE,
                        cex.text = 0.5,
                        cex.lab.x = 0.75,
                        cex.lab.y = 0.65,
                        main=paste("Traits Associated with Metabolomic Modules:", ifelse(sex=="F", "Females", "Males")))
  dev.off()
  
  ### Smaller heatmap with only select traits ###
  
  traitcor.select = traitcor[,c(1,2,4,11:13,15,16)]
  textMatrix.traits.select = textMatrix.traits[,c(1,2,4,11:13,15,16)]
  
  tiff(glue::glue("Trait_Heatmap_selected_traits_{sex}.tiff"), res=400, height=5, width=3.5, units="in")
  map2 = labeledHeatmap(Matrix = traitcor.select,
                        xLabels = colnames(traitcor.select),
                        yLabels = rownames(traitcor.select),
                        ySymbols = rownames(traitcor.select),
                        colorLabels = FALSE,
                        colors = blueWhiteRed(50),
                        textMatrix = textMatrix.traits.select,
                        setStdMargins = FALSE,
                        cex.text = 1,
                        cex.lab.x = 0.75,
                        cex.lab.y = 0.65,
                        main=paste("Traits Associated with Metabolomic Modules:", ifelse(sex=="F", "Females", "Males")))
  dev.off()
  
  
}

lapply(c("F", "M"), traitHeatmaps)



