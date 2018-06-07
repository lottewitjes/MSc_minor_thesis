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

#Change the sample IDs
format_sample_ids = function (atable,method) {
  atable[,1] = as.character(atable[,1])
  for (row in 1:nrow(atable)) {
    sample_id = sprintf("%s_%s", atable[row,1], method)
    atable[row,1] = sample_id
  }
  return(atable)
}

blastx_pfam_count = format_sample_ids(blastx_pfam_count, "blastx")
blastn_pfam_count = format_sample_ids(blastn_pfam_count, "blastn")
blastn_metasapp_pfam_count = format_sample_ids(blastn_metasapp_pfam_count, "metasapp")
blastx_ec_count = format_sample_ids(blastx_ec_count, "blastx")
blastn_ec_count = format_sample_ids(blastn_ec_count, "blastn")
blastn_metasapp_ec_count = format_sample_ids(blastn_metasapp_ec_count, "metasapp")

#Merge the tables
all_pfam_count = rbind(blastx_pfam_count, blastn_pfam_count, blastn_metasapp_pfam_count)
all_ec_count = rbind(blastx_ec_count, blastn_ec_count, blastn_metasapp_ec_count)
  
#Set column names
colnames(all_pfam_count) = c("Sample", "Pfam", "Count")
colnames(all_ec_count) = c("Sample", "EC", "Count")

#Filter for counts
all_pfam_count = all_pfam_count[all_pfam_count$Count > 0,]
all_ec_count = all_ec_count[all_ec_count$Count > 0,]

#Make matrices of the count tables
all_pfam_matrix_count = as.matrix(acast(all_pfam_count, Sample~Pfam, value.var="Count", fill=0))
all_ec_matrix_count = as.matrix(acast(all_ec_count, Sample~EC, value.var="Count", fill=0))

#Filter out Pfam domains/ECs only present in one sample to avoid conflicts with calculating variances
remove_pfam_ec = function(amatrix){
  alist = c()
  for (column in 1:ncol(amatrix)) {
    if (var(amatrix[,column]) != 0 | NA) {
      alist[column] = colnames(amatrix)[column]
    }
  }
  return(alist)
}

#Scale the matrix
all_pfam_matrix_count = scale(all_pfam_matrix_count, center=FALSE, scale=TRUE)
all_ec_matrix_count = scale(all_ec_matrix_count, center=FALSE, scale=TRUE)


method = rep(c("blastn","blastx","metasapp"), 18)
subject = c("1 male","1 male","1 male", 
            "2 female","2 female","2 female", 
            "1 male","1 male","1 male","1 male","1 male","1 male","1 male","1 male", "1 male","1 male","1 male","1 male",
            "2 female","2 female","2 female","2 female","2 female","2 female","2 female","2 female","2 female","2 female","2 female","2 female", 
            "3 male","3 male","3 male","3 male","3 male","3 male","3 male","3 male","3 male","3 male","3 male","3 male", 
            "4 female","4 female","4 female","4 female","4 female","4 female","4 female","4 female","4 female","4 female","4 female","4 female")

#Do the PCA for Pfam domains
all_pfam_pca_count = prcomp(all_pfam_matrix_count, center=FALSE, scale=FALSE)

all_pfam_count_percentage = round(all_pfam_pca_count$sdev / sum(all_pfam_pca_count$sdev) * 100, 2)

all_pfam_count_pc = sprintf("PC%s (%s%%)", which(all_pfam_count_percentage==all_pfam_count_percentage), all_pfam_count_percentage)

all_pfam_pca_count_plot = ggplot(all_pfam_pca_count$x,aes(x=PC1,y=PC2,color=subject)) +
                          geom_point(size=5, aes(shape=method)) +
                          xlab(all_pfam_count_pc[1]) + ylab(all_pfam_count_pc[2]) +
                          labs(title="Pfam domains") +
                          scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                          theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                  axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                  legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                  plot.title=element_text(size=22,colour="black"))
all_pfam_pca_count_plot

#Do the PCA for EC numbers
all_ec_pca_count = prcomp(all_ec_matrix_count, center=FALSE, scale=FALSE)

all_ec_count_percentage = round(all_ec_pca_count$sdev / sum(all_ec_pca_count$sdev) * 100, 2)

all_ec_count_pc = sprintf("PC%s (%s%%)", which(all_ec_count_percentage==all_ec_count_percentage), all_ec_count_percentage)

all_ec_pca_count_plot = ggplot(all_ec_pca_count$x,aes(x=PC1,y=PC2,color=subject)) +
                        geom_point(size=5, aes(shape=method)) +
                        xlab(all_ec_count_pc[1]) + ylab(all_ec_count_pc[2]) +
                        labs(title="ECs") +
                        scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                        theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                plot.title=element_text(size=22,colour="black"))
all_ec_pca_count_plot



