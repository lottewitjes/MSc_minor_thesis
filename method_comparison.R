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
setwd("/media/lottewitjes/Lotte Witjes/MSc_minor_thesis/statistical_analysis/") #Linux
#setwd("D:/MSc_minor_thesis/statistical_analysis/") #Windows

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


Method = rep(c("BLASTN","BLASTX","MetaSAPP"), 18)
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

all_pfam_pca_count_plot = ggplot(all_pfam_pca_count$x,aes(x=PC1,y=PC2,color=Method)) +
                          geom_point(size=5) +
                          xlab(all_pfam_count_pc[1]) + ylab(all_pfam_count_pc[2]) +
                          labs(title="Pfam domains") +
                          scale_color_manual(values=c("#34B233", "#FF7900", "#005172")) +
                          theme_classic() + theme(axis.title.x=element_text(size=24,colour="black"), axis.text.x=element_text(size=22,colour="black"),
                                                  axis.title.y=element_text(size=24,colour="black"), axis.text.y=element_text(size=22,colour="black"),
                                                  legend.title=element_text(size=24,colour="black"), legend.text=element_text(size=22,colour="black"),
                                                  plot.title=element_text(size=24,colour="black"))
all_pfam_pca_count_plot

#Do the PCA for EC numbers
all_ec_pca_count = prcomp(all_ec_matrix_count, center=FALSE, scale=FALSE)

all_ec_count_percentage = round(all_ec_pca_count$sdev / sum(all_ec_pca_count$sdev) * 100, 2)

all_ec_count_pc = sprintf("PC%s (%s%%)", which(all_ec_count_percentage==all_ec_count_percentage), all_ec_count_percentage)

all_ec_pca_count_plot = ggplot(all_ec_pca_count$x,aes(x=PC1,y=PC2,color=Method)) +
                        geom_point(size=5) +
                        xlab(all_ec_count_pc[1]) + ylab(all_ec_count_pc[2]) +
                        labs(title="ECs") +
                        scale_color_manual(values=c("#34B233", "#FF7900", "#005172")) +
                        theme_classic() + theme(axis.title.x=element_text(size=24,colour="black"), axis.text.x=element_text(size=22,colour="black"),
                                                axis.title.y=element_text(size=24,colour="black"), axis.text.y=element_text(size=22,colour="black"),
                                                legend.title=element_text(size=24,colour="black"), legend.text=element_text(size=22,colour="black"),
                                                plot.title=element_text(size=24,colour="black"))
all_ec_pca_count_plot

#Make a Venn-diagram of the three methods
blastx_pfam = unique(all_pfam_count[grep("blastx", all_pfam_count$Sample, fixed=TRUE),2])
blastn_pfam = unique(all_pfam_count[grep("blastn", all_pfam_count$Sample, fixed=TRUE),2])
metasapp_pfam = unique(all_pfam_count[grep("metasapp", all_pfam_count$Sample, fixed=TRUE),2])

blastx_ec = unique(all_ec_count[grep("blastx", all_ec_count$Sample, fixed=TRUE),2])
blastn_ec = unique(all_ec_count[grep("blastn", all_ec_count$Sample, fixed=TRUE),2])
metasapp_ec = unique(all_ec_count[grep("metasapp", all_ec_count$Sample, fixed=TRUE),2])

method_all_pfam = Reduce(intersect, list(blastx_pfam, blastn_pfam, metasapp_pfam)) #core domains
blastx_blastn_pfam = intersect(blastx_pfam, blastn_pfam)
blastx_metasapp_pfam = intersect(blastx_pfam, metasapp_pfam)
blastn_metasapp_pfam = intersect(blastn_pfam, metasapp_pfam)

method_all_ec = Reduce(intersect, list(blastx_ec, blastn_ec, metasapp_ec)) #core domains
blastx_blastn_ec = intersect(blastx_ec, blastn_ec)
blastx_metasapp_ec = intersect(blastx_ec, metasapp_ec)
blastn_metasapp_ec = intersect(blastn_ec, metasapp_ec)

blastx_unique_pfam = Reduce(setdiff, list(blastx_pfam, blastx_blastn_pfam, blastx_metasapp_pfam))
blastn_unique_pfam = Reduce(setdiff, list(blastn_pfam, blastn_metasapp_pfam, blastx_blastn_pfam))
metasapp_unique_pfam = Reduce(setdiff, list(blastn_metasapp_pfam, blastn_metasapp_pfam, blastx_metasapp_pfam))

blastx_unique_ec = Reduce(setdiff, list(blastx_ec, blastx_blastn_ec, blastx_metasapp_ec))
blastn_unique_ec = Reduce(setdiff, list(blastn_ec, blastn_metasapp_ec, blastx_blastn_ec))
metasapp_unique_ec = Reduce(setdiff, list(blastn_metasapp_ec, blastn_metasapp_ec, blastx_metasapp_ec))

grid.newpage()
venn_pfam = draw.triple.venn(area1=length(blastx_pfam), area2=length(blastn_pfam),
                           area3=length(metasapp_pfam),
                           n12=length(blastx_blastn_pfam), n13=length(blastx_metasapp_pfam),
                           n23=length(blastn_metasapp_pfam), n123=length(method_all_pfam),
                           category=c("BLASTX", "BLASTN", "MetaSAPP"),
                           fill=c("#FF7900", "#34B233", "#005172"),
                           cex=2.2, cat.cex=2.2)

grid.newpage()
venn_ec = draw.triple.venn(area1=length(blastx_ec), area2=length(blastn_ec),
                             area3=length(metasapp_ec),
                             n12=length(blastx_blastn_ec), n13=length(blastx_metasapp_ec),
                             n23=length(blastn_metasapp_ec), n123=length(method_all_ec),
                             category=c("BLASTX", "BLASTN", "MetaSAPP"),
                             fill=c("#FF7900", "#34B233", "#005172"),
                             cex=2.2, cat.cex=2.2)



