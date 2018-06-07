#Lotte Witjes
#lottewitjes@outlook.com
#7th of June 2018

library(reshape2)
library(ggplot2)
library(ggfortify)
library(VennDiagram)
library(vegan)
library(dendextend)

#PCA of three methods together###################################################################################################################################
#Set working directory
setwd("/media/lottewitjes/Lotte Witjes/MSc_minor_thesis/statistical_analysis/")

#Load the count tables
blastx_pfam_count = read.table(file="blastx_plots_results/blastx_pfam_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
blastx_ec_count = read.table(file="blastx_plots_results/blastx_ec_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)

blastn_pfam_count = read.table(file="blastn_plots_results/blastn_pfam_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
blastn_ec_count = read.table(file="blastn_plots_results/blastn_ec_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)

blastn_metasapp_pfam_count = read.table(file="blastn_metasapp_plots_results/blastn_metasapp_pfam_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
blastn_metasapp_ec_count = read.table(file="blastn_metasapp_plots_results/blastn_metasapp_ec_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)

#Set column names
colnames(blastx_pfam_count) = c("Sample", "Pfam", "Count")
colnames(blastn_pfam_count) = c("Sample", "Pfam", "Count")
colnames(blastn_metasapp_pfam_count) =c("Sample", "Pfam", "Count")

colnames(blastx_ec_count) = c("Sample", "EC", "Count")
colnames(blastn_ec_count) = c("Sample", "EC", "Count")
colnames(blastn_metasapp_ec_count) = c("Sample", "EC", "Count")

#Filter for counts
blastx_pfam_count = blastx_pfam_count[blastx_pfam_count$Count > 0,]
blastn_pfam_count = blastn_pfam_count[blastn_pfam_count$Count > 0,]
blastn_metasapp_pfam_count = blastn_metasapp_pfam_count[blastn_metasapp_pfam_count$Count > 0,]

blastx_ec_count = blastx_ec_count[blastx_ec_count$Count > 0,]
blastn_ec_count = blastn_ec_count[blastn_ec_count$Count > 0,]
blastn_metasapp_ec_count = blastn_metasapp_ec_count[blastn_metasapp_ec_count$Count > 0,]

format_sample_ids = function (atable) {
  atable[,1] = as.character(atable[,1])
  for (row in 1:nrow(atable)) {
    elements = strsplit(toString(atable[row,1]), "_", fixed=TRUE)
    print(elements)
    elements = elements[[1]]
    sample_id = sprintf("%s_%s", elements[1], elements[2])
    atable[row,1] = sample_id
  }
  return(atable)
}

#Make matrices of the count tables
blastx_pfam_matrix_count = as.matrix(acast(blastx_pfam_count, Sample~Pfam, value.var="Count", fill=0))
blastn_pfam_matrix_count = as.matrix(acast(blastn_pfam_count, Sample~Pfam, value.var="Count", fill=0))
blastn_metasapp_pfam_matrix_count = as.matrix(acast(blastn_metasapp_pfam_count, Sample~Pfam, value.var="Count", fill=0))

blastx_ec_matrix_count = as.matrix(acast(blastx_ec_count, Sample~EC, value.var="Count", fill=0))
blastn_ec_matrix_count = as.matrix(acast(blastn_ec_count, Sample~EC, value.var="Count", fill=0))
blastn_metasapp_ec_matrix_count = as.matrix(acast(blastn_metasapp_ec_count, Sample~EC, value.var="Count", fill=0))

blastx_pfam_matrix_count = cbind(rownames(blastx_pfam_matrix_count), rep("blastx", nrow(blastx_pfam_matrix_count)), blastx_pfam_matrix_count)
blastn_pfam_matrix_count = cbind(rownames(blastn_pfam_matrix_count), rep("blastn", nrow(blastn_pfam_matrix_count)), blastn_pfam_matrix_count)
blastn_metasapp_pfam_matrix_count = cbind(rownames(blastn_metasapp_pfam_matrix_count), rep("blastn_metasapp", nrow(blastn_metasapp_pfam_matrix_count)), blastn_metasapp_pfam_matrix_count)
blastx_ec_matrix_count = cbind(rownames(blastx_ec_matrix_count), rep("blastx", nrow(blastx_ec_matrix_count)), blastx_ec_matrix_count)
blastn_ec_matrix_count = cbind(rownames(blastn_ec_matrix_count), rep("blastn", nrow(blastn_ec_matrix_count)), blastn_ec_matrix_count)
blastn_metasapp_ec_matrix_count = cbind(rownames(blastn_metasapp_ec_matrix_count), rep("blastn_metasapp", nrow(blastn_metasapp_ec_matrix_count)), blastn_metasapp_ec_matrix_count)

all_pfam_matrix_count = as.matrix(merge(blastx_pfam_matrix_count, blastn_pfam_matrix_count, all=TRUE))
all_pfam_matrix_count = as.matrix(merge(all_pfam_matrix_count, blastn_metasapp_pfam_matrix_count, all=TRUE))
rownames(all_pfam_matrix_count) = all_pfam_matrix_count[,1]
all_pfam_matrix_count = all_pfam_matrix_count[,-1]
method_names = all_pfam_matrix_count[,1]
sample_names = rownames(all_pfam_matrix_count)
all_pfam_matrix_count = all_pfam_matrix_count[,-1]
storage.mode(all_pfam_matrix_count) = "numeric"
all_pfam_matrix_count = all_pfam_matrix_count[,apply(all_pfam_matrix_count, 2, var) != NA]

all_pfam_pca_count = prcomp(na.omit(all_pfam_matrix_count), center=FALSE, scale=FALSE)

all_pfam_count_percentage = round(all_pfam_pca_count$sdev / sum(all_pfam_pca_count$sdev) * 100, 2)

all_pfam_count_pc = sprintf("PC%s (%s%%)", which(all_pfam_count_percentage==all_pfam_count_percentage), all_pfam_count_percentage)

all_pfam_pca_count_plot = ggplot(all_pfam_pca_count$x,aes(x=PC1,y=PC2,color=rownames(all_pfam_matrix_count))) +
  geom_point(size=5) +
  xlab(all_pfam_count_pc[1]) + ylab(all_pfam_count_pc[2]) +
  labs(title="Pfam domains") +
  scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
  theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                          axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                          legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                          plot.title=element_text(size=22,colour="black"))
all_pfam_pca_count_plot
