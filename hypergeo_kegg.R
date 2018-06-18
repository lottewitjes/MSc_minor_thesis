#Lotte Witjes
#lottewitjes@outlook.com
#7th of June 2018

library(reshape2)
library(ggplot2)
library(ggfortify)
library(VennDiagram)
library(vegan)
library(dendextend)

#Hypergeometric test on KEGG pathways###########################################################################################################################
setwd("/media/lottewitjes/Lotte Witjes/MSc_minor_thesis/statistical_analysis/") #Linux
setwd("D:/MSc_minor_thesis/statistical_analysis/") #Windows

blastn_kegg_counts = read.table(file="blastn_plots_results/blastn_kegg_enzyme_count.tsv", sep="\t", header=TRUE)
blastn_kegg_counts = cbind(blastn_kegg_counts, rep(794, length(nrow(blastn_kegg_counts))), rep(3440, length(nrow(blastn_kegg_counts))))
colnames(blastn_kegg_counts) = c("pathway", "mapped_enzymes", "enzymes_in_pathway", "total_sample", "total_ec")

hypgeo_kegg = function(counts) {
  p_values = c()
  for (row in 1:nrow(counts)) {
    if (counts[row,2] > 0) {
      hypgeo = 1-phyper(counts[row,2], counts[row,3], (counts[row,5]-counts[row,3]), counts[row,4])
      p_values[row] = hypgeo
    } else {p_values[row] = NA}
  }
  return(p_values)
}

blastn_kegg_counts = cbind(blastn_kegg_counts, hypgeo_kegg(blastn_kegg_counts))
colnames(blastn_kegg_counts)[6] = "p_value_hypgeo"
blastn_kegg_counts$p_value_hypgeo = p.adjust(blastn_kegg_counts$p_value_hypgeo, method="BH")

blastn_kegg_counts = blastn_kegg_counts[order(blastn_kegg_counts$p_value_hypgeo),]
blastn_kegg_counts = blastn_kegg_counts[which(blastn_kegg_counts$p_value_hypgeo < 0.05),]
