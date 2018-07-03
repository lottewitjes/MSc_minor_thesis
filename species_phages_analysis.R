#Lotte Witjes
#lottewitjes@outlook.com
#3rd of July 2018

library(reshape2)
library(ggplot2)
library(ggfortify)
library(VennDiagram)
library(vegan)
library(dendextend)

#Data preprocessing#####################################################################################################################################################
#Set working directory
setwd("/media/lottewitjes/Lotte Witjes/MSc_minor_thesis/statistical_analysis/blastn_results/")

#Load the data
species_counts = read.table(file="blastn_species_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
phages_counts = read.table(file="blastn_phage_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)

#Set column names
colnames(species_counts) = c("Sample", "Species", "Count")
colnames(phages_counts) = c("Sample", "Phage", "Count")

#Filter on counts
species_counts = species_counts[species_counts$Count > 5,]
phages_counts = phages_counts[phages_counts$Count > 5,]

#Specify the labels
sample_names = c("m.s1.d3.mrn.rep", "f.s2.d1.mrn.rep", "m.s1.d1.mrn", "m.s1.d1.aft", "m.s1.d3.mrn", "m.s1.d3.aft",
                 "f.s2.d1.mrn", "f.s2.d1.aft", "f.s2.d3.mrn", "f.s2.d3.aft", "m.s3.d1.mrn", "m.s3.d1.aft", "m.s3.d3.mrn",
                 "m.s3.d3.aft", "f.s4.d1.mrn", "f.s4.d1.aft", "f.s4.d3.mrn", "f.s4.d3.aft")
subject = c("1 male", "2 female", "1 male", "1 male", "1 male", "1 male", "2 female", "2 female", "2 female", "2 female", "3 male", "3 male",
            "3 male", "3 male", "4 female", "4 female", "4 female", "4 female")

#Transform data into matrices
species_matrix_count = as.matrix(acast(species_counts, Sample~Species, value.var="Count", fill=0))
phages_matrix_count = as.matrix(acast(phages_counts, Sample~Phage, value.var="Count", fill=0))

#Scale the matrices
species_matrix_count_scaled = scale(species_matrix_count, center=FALSE, scale=TRUE)
phages_matrix_count_scaled = scale(phages_matrix_count, center=FALSE, scale=TRUE)

#Beta diversity and richness plots for species and phages###############################################################################################################
beta_diversity_species = diversity(species_matrix_count_scaled, index="shannon")
richness_species = specnumber(species_matrix_count)

beta_diversity_phages = diversity(phages_matrix_count_scaled, index="shannon")
richness_phages = specnumber(phages_matrix_count)

timepoints = c("Day 3, morning", "Day 1, morning", "Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon",
               "Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon", "Day 1, morning", "Day 1, afternoon",
               "Day 3, morning", "Day 3, afternoon", "Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")

beta_diversity_df = data.frame(subject, timepoints, beta_diversity_species, beta_diversity_phages, richness_species, richness_phages)
beta_diversity_df = beta_diversity_df[-c(1,2),]

beta_diversity_species_plot = ggplot(data=beta_diversity_df, aes(x=timepoints, y=beta_diversity_species, group=subject)) +
                              geom_line(aes(color=subject), size=3) +
                              geom_point(aes(color=subject), size=5) +
                              scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                              scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                              xlab("Timepoints") + ylab("Beta diversity (Shannon Index)") +
                              labs(title="Species") +
                              theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                      axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                      legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                      plot.title=element_text(size=22,colour="black"))
beta_diversity_species_plot

beta_diversity_phages_plot = ggplot(data = beta_diversity_df, aes(x=timepoints, y=beta_diversity_phages, group=subject)) +
                             geom_line(aes(color=subject), size=3) +
                             geom_point(aes(color=subject), size=5) +
                             scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                             scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                             xlab("Timepoints") + ylab("Beta diversity (Shannon Index)") +
                             labs(title="Phages") +
                             theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                     axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                     legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                     plot.title=element_text(size=22,colour="black"))
beta_diversity_phages_plot

richness_species_plot = ggplot(data = beta_diversity_df, aes(x=timepoints, y=richness_species, group=subject)) +
                        geom_line(aes(color=subject), size=3) +
                        geom_point(aes(color=subject), size=5) +
                        scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                        scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                        xlab("Timepoints") + ylab("Richness") +
                        labs(title="Species") +
                        theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                plot.title=element_text(size=22,colour="black"))
richness_species_plot

richness_phages_plot = ggplot(data = beta_diversity_df, aes(x=timepoints, y=richness_phages, group=subject)) +
                       geom_line(aes(color=subject), size=3) +
                       geom_point(aes(color=subject), size=5) +
                       scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                       scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                       xlab("Timepoints") + ylab("Richness") +
                       labs(title="Phages") +
                       theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                               axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                               legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                               plot.title=element_text(size=22,colour="black"))
richness_phages_plot

#Bar chart representing genus for each sample#########################################################################################################################
average_species_counts = colMeans(species_matrix_count)
average_species_counts = sort(average_species_counts, decreasing=TRUE)
top10_species = names(average_species_counts[1:10])

`%ni%` <- Negate(`%in%`)

other_species = species_counts[species_counts$Species %ni% top10_species,]
other_sum = aggregate(other_species$Count, by=list(Category=other_species$Sample), FUN=sum)
other_sum = cbind(other_sum[,1], rep("Other", 18), other_sum[,2])
colnames(other_sum) = c("Sample", "Species", "Count")

top10_species = c(top10_species, "Other")
species_counts_other = rbind(species_counts, other_sum)
species_counts_subset = species_counts_other[species_counts_other$Species %in% top10_species,]

part_day = c("Mrn", "Mrn", 
             "Mrn", "Aft", "Mrn", "Aft",
             "Mrn", "Aft", "Mrn", "Aft",
             "Mrn", "Aft", "Mrn", "Aft",
             "Mrn", "Aft", "Mrn", "Aft")

which_day = c("Day 3", "Day 1", 
              "Day 1", "Day 1", "Day 3", "Day 3",
              "Day 1", "Day 1", "Day 3", "Day 3", 
              "Day 1", "Day 1","Day 3", "Day 3", 
              "Day 1", "Day 1", "Day 3", "Day 3")
subject = c("1 male", "2 female", 
            "1 male", "1 male", "1 male", "1 male", 
            "2 female", "2 female", "2 female", "2 female", 
            "3 male", "3 male","3 male", "3 male", 
            "4 female", "4 female", "4 female", "4 female")


bar_plot_species = ggplot(data=species_counts_subset, aes(fill=Species, y=as.numeric(Count), x=Sample)) +
                   geom_bar(stat="identity", position="fill") +
                   scale_x_discrete(labels=sample_names) + coord_flip() +
                   xlab("Samples") + ylab("Relative contribution") +
                   labs(title="Specific active species (mRNA)") +
                   scale_color_brewer(type="qual", palette = 2) +
                   theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=13,colour="black"),
                                           axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                           legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                           legend.position="bottom", plot.title=element_text(size=22,colour="black"))
bar_plot_species

#Top X most present phages and species for each sample################################################################################################################
top_x_phages = data.frame(row.names=c("1", "2", "3", "4", "5", "6", "7", "8", "9","10"))
for (row in 1:nrow(phages_matrix_count)) {
  top_x = sort(phages_matrix_count[row,], decreasing=TRUE)
  top_x_phages = cbind(top_x_phages, names(top_x[1:10]))
}
colnames(top_x_phages) = sample_names

top_x_species = data.frame(row.names=c("1", "2", "3", "4", "5", "6", "7", "8", "9","10"))
for (row in 1:nrow(species_matrix_count)) {
  top_x = sort(species_matrix_count[row,], decreasing=TRUE)
  top_x_species = cbind(top_x_species, names(top_x[1:10]))
}
colnames(top_x_species) = sample_names
