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
setwd("D:/MSc_minor_thesis/statistical_analysis/species_genera_phages/")

#Load the data
species_counts = read.table(file="blastn_species_counts_without_strain.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
phages_counts = read.table(file="blastn_phage_count.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)

#Set column names
colnames(species_counts) = c("Sample", "Species", "Count")
colnames(phages_counts) = c("Sample", "Phage", "Count")

#Filter on counts
species_counts = species_counts[species_counts$Count > 0,]
phages_counts = phages_counts[phages_counts$Count > 0,]

#Specify the labels
sample_names = c("m.s1.d1.mrn", "m.s1.d1.aft", "m.s1.d3.mrn", "m.s1.d3.aft",
                 "f.s2.d1.mrn", "f.s2.d1.aft", "f.s2.d3.mrn", "f.s2.d3.aft", 
                 "m.s3.d1.mrn", "m.s3.d1.aft", "m.s3.d3.mrn", "m.s3.d3.aft", 
                 "f.s4.d1.mrn", "f.s4.d1.aft", "f.s4.d3.mrn", "f.s4.d3.aft")
subject = c("1 male", "1 male", "1 male", "1 male", 
            "2 female", "2 female", "2 female", "2 female", 
            "3 male", "3 male", "3 male", "3 male", 
            "4 female", "4 female", "4 female", "4 female")

#Transform data into matrices
species_matrix_count = as.matrix(acast(species_counts, Sample~Species, value.var="Count", fill=0, fun.aggregate=sum))
phages_matrix_count = as.matrix(acast(phages_counts, Sample~Phage, value.var="Count", fill=0, fun.aggregate=sum))

#Scale the matrices
species_matrix_count_scaled = scale(species_matrix_count, center=FALSE, scale=TRUE)
phages_matrix_count_scaled = scale(phages_matrix_count, center=FALSE, scale=TRUE)

#Shannon diversity and richness plots for species and phages###############################################################################################################
shannon_diversity_species = diversity(species_matrix_count, index="shannon")
richness_species = specnumber(species_matrix_count)

shannon_diversity_phages = diversity(phages_matrix_count, index="shannon")
richness_phages = specnumber(phages_matrix_count)

timepoints = c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon",
               "Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon", "Day 1, morning", "Day 1, afternoon",
               "Day 3, morning", "Day 3, afternoon", "Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")

shannon_diversity_df = data.frame(subject, timepoints, shannon_diversity_species, shannon_diversity_phages, richness_species, richness_phages)

shannon_diversity_species_plot = ggplot(data=shannon_diversity_df, aes(x=timepoints, y=shannon_diversity_species, group=subject)) +
                              geom_line(aes(color=subject), size=3) +
                              geom_point(aes(color=subject), size=5) +
                              scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                              scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                              xlab("Timepoints") + ylab("Shannon diversity") +
                              labs(title="Species") +
                              theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                      axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                      legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                      plot.title=element_text(size=22,colour="black"))
shannon_diversity_species_plot

shannon_diversity_phages_plot = ggplot(data = shannon_diversity_df, aes(x=timepoints, y=shannon_diversity_phages, group=subject)) +
                             geom_line(aes(color=subject), size=3) +
                             geom_point(aes(color=subject), size=5) +
                             scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                             scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                             xlab("Timepoints") + ylab("Shannon diversity") +
                             labs(title="Phages") +
                             theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                     axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                     legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                     plot.title=element_text(size=22,colour="black"))
shannon_diversity_phages_plot

richness_species_plot = ggplot(data = shannon_diversity_df, aes(x=timepoints, y=richness_species, group=subject)) +
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

richness_phages_plot = ggplot(data = shannon_diversity_df, aes(x=timepoints, y=richness_phages, group=subject)) +
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
genus_counts = read.table(file="blastn_genus_counts.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(genus_counts) = c("Sample", "Genus", "Count")
genus_counts = genus_counts[genus_counts$Count > 0,]
genus_matrix_count = as.matrix(acast(genus_counts, Sample~Genus, value.var="Count", fill=0, fun.aggregate=sum))

average_genus_counts = colMeans(genus_matrix_count)
average_genus_counts = sort(average_genus_counts, decreasing=TRUE)
top10_genus = names(average_genus_counts[1:9])

`%ni%` <- Negate(`%in%`)

other_genus = genus_counts[genus_counts$Genus %ni% top10_genus,]
other_sum = aggregate(other_genus$Count, by=list(Category=other_genus$Sample), FUN=sum)
other_sum = cbind(other_sum[,1], rep("Other", 16), other_sum[,2])
colnames(other_sum) = c("Sample", "Genus", "Count")

top10_genus = c(top10_genus, "Other")
genus_counts_other = rbind(genus_counts, other_sum)
genus_counts_subset = genus_counts_other[genus_counts_other$Genus %in% top10_genus,]

part_day = c("Mrn", "Aft", "Mrn", "Aft",
             "Mrn", "Aft", "Mrn", "Aft",
             "Mrn", "Aft", "Mrn", "Aft",
             "Mrn", "Aft", "Mrn", "Aft")

which_day = c("Day 1", "Day 1", "Day 3", "Day 3",
              "Day 1", "Day 1", "Day 3", "Day 3", 
              "Day 1", "Day 1","Day 3", "Day 3", 
              "Day 1", "Day 1", "Day 3", "Day 3")

subject = c("1 male", "1 male", "1 male", "1 male", 
            "2 female", "2 female", "2 female", "2 female", 
            "3 male", "3 male","3 male", "3 male", 
            "4 female", "4 female", "4 female", "4 female")


bar_plot_genus = ggplot(data=genus_counts_subset, aes(fill=Genus, y=as.numeric(Count), x=Sample)) +
                 geom_bar(stat="identity", position="fill") +
                 scale_x_discrete(labels=sample_names) + coord_flip() +
                 scale_y_continuous(expand=c(0,0)) +
                 xlab("Samples") + ylab("Relative contribution") +
                 labs(title="Active fraction (16S rRNA)") +
                 scale_fill_brewer(palette="Spectral") +
                 theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=13,colour="black"),
                                         axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                         legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                         legend.position="bottom", plot.title=element_text(size=22,colour="black"))
bar_plot_genus

#PCA with count data############################################################################################################################################
species_pca_count = prcomp(species_matrix_count_scaled, center=FALSE, scale=FALSE)
phages_pca_count = prcomp(phages_matrix_count_scaled, center=FALSE, scale=FALSE)

species_count_percentage = round(species_pca_count$sdev / sum(species_pca_count$sdev) * 100, 2)
phages_count_percentage = round(phages_pca_count$sdev / sum(phages_pca_count$sdev) * 100, 2)

species_count_pc = sprintf("PC%s (%s%%)", which(species_count_percentage==species_count_percentage), species_count_percentage)
phages_count_pc = sprintf("PC%s (%s%%)", which(phages_count_percentage==phages_count_percentage), phages_count_percentage)

species_pca_count_plot = ggplot(species_pca_count$x,aes(x=PC1,y=PC2,color=subject)) +
                         geom_point(size=5) +
                         xlab(species_count_pc[1]) + ylab(species_count_pc[2]) +
                         labs(title="Species") +
                         scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                         theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                 axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                 legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                 plot.title=element_text(size=22,colour="black"))
species_pca_count_plot

phages_pca_count_plot = ggplot(phages_pca_count$x,aes(x=PC1,y=PC2,color=subject)) +
                        geom_point(size=5) +
                        xlab(phages_count_pc[1]) + ylab(phages_count_pc[2]) +
                        labs(title="Phages") +
                        scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                        theme_classic() + theme(axis.title.x=element_text(size=22,colour="black"), axis.text.x=element_text(size=20,colour="black"),
                                                axis.title.y=element_text(size=22,colour="black"), axis.text.y=element_text(size=20,colour="black"),
                                                legend.title=element_text(size=22,colour="black"), legend.text=element_text(size=20,colour="black"),
                                                plot.title=element_text(size=22,colour="black"))
phages_pca_count_plot

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

#Venn-diagram BETWEEN subjects##################################################################################################################################
#Obtain list of species in each subject
subject1_species = unique(species_counts[species_counts$Sample %in% c("NG-5450_A", "NG-5593_1A", "NG-5593_1B", "NG-5593_1C", "NG-5593_1D"),2])
subject2_species = unique(species_counts[species_counts$Sample %in% c("NG-5450_B", "NG-5593_2A", "NG-5593_2B", "NG-5593_2C", "NG-5593_2D"),2])
subject3_species = unique(species_counts[species_counts$Sample %in% c("NG-5593_3A", "NG-5593_3B", "NG-5593_3C", "NG-5593_3D"),2])
subject4_species = unique(species_counts[species_counts$Sample %in% c("NG-5593_4A", "NG-5593_4B", "NG-5593_4C", "NG-5593_4D"),2])

subject1_phages = unique(phages_counts[phages_counts$Sample %in% c("NG-5450_A", "NG-5593_1A", "NG-5593_1B", "NG-5593_1C", "NG-5593_1D"),2])
subject2_phages = unique(phages_counts[phages_counts$Sample %in% c("NG-5450_B", "NG-5593_2A", "NG-5593_2B", "NG-5593_2C", "NG-5593_2D"),2])
subject3_phages = unique(phages_counts[phages_counts$Sample %in% c("NG-5593_3A", "NG-5593_3B", "NG-5593_3C", "NG-5593_3D"),2])
subject4_phages = unique(phages_counts[phages_counts$Sample %in% c("NG-5593_4A", "NG-5593_4B", "NG-5593_4C", "NG-5593_4D"),2])

#Find species and phages present or absent in all and between subjects
subject_all_species = Reduce(intersect, list(subject1_species, subject2_species, subject3_species, subject4_species)) #core species
subject_1_2_species = intersect(subject1_species, subject2_species)
subject_1_3_species = intersect(subject1_species, subject3_species)
subject_1_4_species = intersect(subject1_species, subject4_species)
subject_2_3_species = intersect(subject2_species, subject3_species)
subject_2_4_species = intersect(subject2_species, subject4_species)
subject_3_4_species = intersect(subject3_species, subject4_species)
subject_1_2_3_species = Reduce(intersect, list(subject1_species, subject2_species, subject3_species))
subject_1_2_4_species = Reduce(intersect, list(subject1_species, subject2_species, subject4_species))
subject_1_3_4_species = Reduce(intersect, list(subject1_species, subject3_species, subject4_species))
subject_2_3_4_species = Reduce(intersect, list(subject2_species, subject3_species, subject4_species))

subject_all_phages = Reduce(intersect, list(subject1_phages, subject2_phages, subject3_phages, subject4_phages)) #core phages
subject_1_2_phages = intersect(subject1_phages, subject2_phages)
subject_1_3_phages = intersect(subject1_phages, subject3_phages)
subject_1_4_phages = intersect(subject1_phages, subject4_phages)
subject_2_3_phages = intersect(subject2_phages, subject3_phages)
subject_2_4_phages = intersect(subject2_phages, subject4_phages)
subject_3_4_phages = intersect(subject3_phages, subject4_phages)
subject_1_2_3_phages = Reduce(intersect, list(subject1_phages, subject2_phages, subject3_phages))
subject_1_2_4_phages = Reduce(intersect, list(subject1_phages, subject2_phages, subject4_phages))
subject_1_3_4_phages = Reduce(intersect, list(subject1_phages, subject3_phages, subject4_phages))
subject_2_3_4_phages = Reduce(intersect, list(subject2_phages, subject3_phages, subject4_phages))

#Venn-diagram of subjects
grid.newpage()
venn_species = draw.quad.venn(area1=length(subject1_species), area2=length(subject2_species),
                           area3=length(subject3_species), area4=length(subject4_species),
                           n12=length(subject_1_2_species), n13=length(subject_1_3_species),
                           n14=length(subject_1_4_species), n23=length(subject_2_3_species),
                           n24=length(subject_2_4_species), n34=length(subject_3_4_species),
                           n1234=length(subject_all_species), n123=length(subject_1_2_3_species),
                           n124=length(subject_1_2_4_species), n134=length(subject_1_3_4_species),
                           n234=length(subject_2_3_4_species),
                           category=c("Subject 1 male", "Subject 2 female", 
                                      "Subject 3 male", "Subject 4 female"),
                           fill=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c"),
                           cex=2.2, cat.cex=2.2)

grid.newpage()
venn_phages = draw.quad.venn(area1=length(subject1_phages), area2=length(subject2_phages),
                           area3=length(subject3_phages), area4=length(subject4_phages),
                           n12=length(subject_1_2_phages), n13=length(subject_1_3_phages),
                           n14=length(subject_1_4_phages), n23=length(subject_2_3_phages),
                           n24=length(subject_2_4_phages), n34=length(subject_3_4_phages),
                           n1234=length(subject_all_phages), n123=length(subject_1_2_3_phages),
                           n124=length(subject_1_2_4_phages), n134=length(subject_1_3_4_phages),
                           n234=length(subject_2_3_4_phages),
                           category=c("Subject 1 male", "Subject 2 female", 
                                      "Subject 3 male", "Subject 4 female"),
                           fill=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c"),
                           cex=2.2, cat.cex=2.2)

#Unique species and phages per subject####################################################################################################
subject1_unique_species = Reduce(setdiff, list(subject1_species, subject_1_2_species, subject_1_2_3_species, subject_1_3_species, subject_1_4_species,
                                            subject_1_2_4_species, subject_all_species, subject_1_3_4_species))
subject2_unique_species = Reduce(setdiff, list(subject2_species, subject_1_2_species, subject_1_2_3_species, subject_2_3_species, subject_2_4_species,
                                            subject_1_2_4_species, subject_all_species, subject_2_3_4_species))
subject3_unique_species = Reduce(setdiff, list(subject3_species, subject_1_3_species, subject_1_2_3_species, subject_2_3_species, subject_3_4_species,
                                            subject_all_species, subject_2_3_4_species))
subject4_unique_species = Reduce(setdiff, list(subject4_species, subject_1_4_species, subject_2_4_species, subject_3_4_species, subject_1_3_4_species,
                                            subject_1_2_4_species, subject_2_3_4_species))

subject1_unique_phages = Reduce(setdiff, list(subject1_phages, subject_1_2_phages, subject_1_2_3_phages, subject_1_3_phages, subject_1_4_phages,
                                          subject_1_2_4_phages, subject_all_phages, subject_1_3_4_phages))
subject2_unique_phages = Reduce(setdiff, list(subject2_phages, subject_1_2_phages, subject_1_2_3_phages, subject_2_3_phages, subject_2_4_phages,
                                          subject_1_2_4_phages, subject_all_phages, subject_2_3_4_phages))
subject3_unique_phages = Reduce(setdiff, list(subject3_phages, subject_1_3_phages, subject_1_2_3_phages, subject_2_3_phages, subject_3_4_phages,
                                          subject_all_phages, subject_2_3_4_phages))
subject4_unique_phages = Reduce(setdiff, list(subject4_phages, subject_1_4_phages, subject_2_4_phages, subject_3_4_phages, subject_1_3_4_phages,
                                          subject_1_2_4_phages, subject_2_3_4_phages))

subject1_samples = c("NG-5593_1A", "NG-5593_1B", "NG-5593_1C", "NG-5593_1D")
subject2_samples = c("NG-5593_2A", "NG-5593_2B", "NG-5593_2C", "NG-5593_2D")
subject3_samples = c("NG-5593_3A", "NG-5593_3B", "NG-5593_3C", "NG-5593_3D")
subject4_samples = c("NG-5593_4A", "NG-5593_4B", "NG-5593_4C", "NG-5593_4D")

unique_species_sub1 = species_matrix_count[subject1_samples,subject1_unique_species]
unique_species_sub1 = colMeans(unique_species_sub1)
unique_species_sub1 = sort(unique_species_sub1, decreasing=TRUE)

unique_species_sub2 = species_matrix_count[subject2_samples,subject2_unique_species]
unique_species_sub2 = colMeans(unique_species_sub2)
unique_species_sub2 = sort(unique_species_sub2, decreasing=TRUE)

unique_species_sub3 = species_matrix_count[subject3_samples,subject3_unique_species]
unique_species_sub3 = colMeans(unique_species_sub3)
unique_species_sub3 = sort(unique_species_sub3, decreasing=TRUE)

unique_species_sub4 = species_matrix_count[subject4_samples,subject4_unique_species]
unique_species_sub4 = colMeans(unique_species_sub4)
unique_species_sub4 = sort(unique_species_sub4, decreasing=TRUE)

unique_phages_sub1 = phages_matrix_count[subject1_samples,subject1_unique_phages]
unique_phages_sub1 = colMeans(unique_phages_sub1)
unique_phages_sub1 = sort(unique_phages_sub1, decreasing=TRUE)

unique_phages_sub2 = phages_matrix_count[subject2_samples,subject2_unique_phages]
unique_phages_sub2 = colMeans(unique_phages_sub2)
unique_phages_sub2 = sort(unique_phages_sub2, decreasing=TRUE)

unique_phages_sub3 = phages_matrix_count[subject3_samples,subject3_unique_phages]
unique_phages_sub3 = colMeans(unique_phages_sub3)
unique_phages_sub3 = sort(unique_phages_sub3, decreasing=TRUE)

unique_phages_sub4 = phages_matrix_count[subject4_samples,subject4_unique_phages]
unique_phages_sub4 = colMeans(unique_phages_sub4)
unique_phages_sub4 = sort(unique_phages_sub4, decreasing=TRUE)

core_species = species_matrix_count[, subject_all_species]
core_species = colMeans(core_species)
core_species = sort(core_species, decreasing=TRUE)

core_phages = phages_matrix_count[, subject_all_phages]
core_phages = colMeans(core_phages)
core_phages = sort(core_phages, decreasing=TRUE)