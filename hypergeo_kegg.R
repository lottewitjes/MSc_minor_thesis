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
setwd("/media/lottewitjes/Lotte Witjes/MSc_minor_thesis/statistical_analysis/")

blastx_kegg_counts = read.table(file="blastx_plots_results/blastx_kegg_enzyme_count.tsv", sep="\t", header=TRUE)
blastx_kegg_counts = cbind(blastx_kegg_counts, rep(641, length(nrow(blastx_kegg_counts))), rep(ncol(ec_matrix_count), length(nrow(blastx_kegg_counts))))
colnames(blastx_kegg_counts) = c("pathway", "mapped_enzymes", "enzymes_in_pathway", "total_core", "total_ec")

blastn_kegg_counts = read.table(file="blastn_plots_results/blastn_kegg_enzyme_count.tsv", sep="\t", header=TRUE)
blastn_kegg_counts = cbind(blastn_kegg_counts, rep(924, length(nrow(blastn_kegg_counts))), rep(ncol(ec_matrix_count), length(nrow(blastn_kegg_counts))))
colnames(blastn_kegg_counts) = colnames(blastx_kegg_counts)

blastn_metasapp_kegg_counts = read.table(file="blastn_metasapp_plots_results/blastn_metasapp_kegg_enzyme_count.tsv", sep="\t", header=TRUE)
blastn_metasapp_kegg_counts = cbind(blastn_metasapp_kegg_counts, rep(1164, length(nrow(blastn_metasapp_kegg_counts))), rep(ncol(ec_matrix_count), length(nrow(blastn_metasapp_kegg_counts))))
colnames(blastn_metasapp_kegg_counts) = colnames(blastx_kegg_counts)

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

blastx_kegg_counts = cbind(blastx_kegg_counts, hypgeo_kegg(blastx_kegg_counts))
colnames(blastx_kegg_counts)[6] = "p_value_hypgeo"
blastx_kegg_counts = blastx_kegg_counts[order(blastx_kegg_counts$p_value_hypgeo),]
blastx_kegg_counts = blastx_kegg_counts[which(blastx_kegg_counts$p_value_hypgeo < 0.05),]

blastn_kegg_counts = cbind(blastn_kegg_counts, hypgeo_kegg(blastn_kegg_counts))
colnames(blastn_kegg_counts)[6] = "p_value_hypgeo"
blastn_kegg_counts = blastn_kegg_counts[order(blastn_kegg_counts$p_value_hypgeo),]
blastn_kegg_counts = blastn_kegg_counts[which(blastn_kegg_counts$p_value_hypgeo < 0.05),]

blastn_metasapp_kegg_counts = cbind(blastn_metasapp_kegg_counts, hypgeo_kegg(blastn_metasapp_kegg_counts))
colnames(blastn_metasapp_kegg_counts)[6] = "p_value_hypgeo"
blastn_metasapp_kegg_counts = blastn_metasapp_kegg_counts[order(blastn_metasapp_kegg_counts$p_value_hypgeo),]
blastn_metasapp_kegg_counts = blastn_metasapp_kegg_counts[which(blastn_metasapp_kegg_counts$p_value_hypgeo < 0.05),]