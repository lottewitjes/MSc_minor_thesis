#Lotte Witjes
#lottewitjes@outlook.com
#22 of May 2018

#Load libraries
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("ggfortify")
#install.packages("VennDiagram")
#install.packages("vegan")
#install.packages("dendextend")
library(reshape2)
library(ggplot2)
library(ggfortify)
library(VennDiagram)
library(vegan)
library(dendextend)

#Set working directory
setwd("/media/lottewitjes/Lotte Witjes/MSc_minor_thesis/statistical_analysis/")

##################################reads#########################################################################################################################
#Load data
read_data = read.table(file="data_overview.tsv", sep="\t", header=TRUE)
read_data = read_data[,1:5]
read_data_tp = melt(read_data)
read_data_tp$X = factor(read_data_tp$X, levels=c("m.s1.d1.mrn", "m.s1.d1.aft", "m.s1.d3.mrn", "m.s1.d3.mrn.rep", "m.s1.d3.aft",
                                                 "f.s2.d1.mrn", "f.s2.d1.mrn.rep", "f.s2.d1.aft", "f.s2.d3.mrn", "f.s2.d3.aft", "m.s3.d1.mrn", "m.s3.d1.aft", "m.s3.d3.mrn",
                                                 "m.s3.d3.aft", "f.s4.d1.mrn", "f.s4.d1.aft", "f.s4.d3.mrn", "f.s4.d3.aft"))
read_data_tp$variable = factor(read_data_tp$variable, levels=c("Nonhuman.mRNA", "Human.mRNA", "rRNA", "Adaptor.and.low.quality"))

read_plot = ggplot(data=read_data_tp, aes(x=X, y=value, fill=variable)) +
            geom_bar(stat="identity") +
            xlab("Sample") + ylab("Number of reads") + coord_flip() +
            labs(title="Reads overview per sample") +
            scale_y_continuous(expand=c(0,0),labels=scales::comma) +
            scale_fill_brewer(palette="Paired", name="Type of reads", breaks=c("Nonhuman.mRNA", "Human.mRNA", "rRNA", "Adaptor.and.low.quality"),
                              labels=c("Non-human mRNA", "Human mRNA", "rRNA", "Adaptor and low quality")) +
            theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                    axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                    legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                    plot.title=element_text(size=20,colour="black"))
read_plot

###############################domain/enzyme analyses and plots#################################################################################################
#Load data
setwd("/media/lottewitjes/Lotte Witjes/MSc_minor_thesis/statistical_analysis/blastx_plots_results/")
pfam_count = read.table(file="blastx_pfam_count.tsv", sep="\t", header=FALSE)
ec_count = read.table(file="blastx_ec_count.tsv", sep="\t", header=FALSE)

#Set column names
colnames(pfam_count) = c("Sample", "Pfam", "Count")
colnames(ec_count) = c("Sample", "EC", "Count")

#Filter on counts
pfam_count = pfam_count[pfam_count$Count > 5,]
ec_count = ec_count[ec_count$Count > 5,]

#Specify the labels
sample_names = c("m.s1.d3.mrn.rep", "f.s2.d1.mrn.rep", "m.s1.d1.mrn", "m.s1.d1.aft", "m.s1.d3.mrn", "m.s1.d3.aft",
                       "f.s2.d1.mrn", "f.s2.d1.aft", "f.s2.d3.mrn", "f.s2.d3.aft", "m.s3.d1.mrn", "m.s3.d1.aft", "m.s3.d3.mrn",
                       "m.s3.d3.aft", "f.s4.d1.mrn", "f.s4.d1.aft", "f.s4.d3.mrn", "f.s4.d3.aft")
subject = c("1 male", "2 female", "1 male", "1 male", "1 male", "1 male", "2 female", "2 female", "2 female", "2 female", "3 male", "3 male",
            "3 male", "3 male", "4 female", "4 female", "4 female", "4 female")

#Transform data into matrices
pfam_matrix_count = as.matrix(acast(pfam_count, Sample~Pfam, value.var="Count", fill=0))

ec_matrix_count = as.matrix(acast(ec_count, Sample~EC, value.var="Count", fill=0))

#Domain/enzyme alpha/beta/gamma diversity#######################################################################################################################
beta_diversity_pfam = diversity(pfam_matrix_count, index="shannon")
alpha_diversity_pfam = fisher.alpha(pfam_matrix_count)
richness_pfam = specnumber(pfam_matrix_count)

beta_diversity_ec = diversity(ec_matrix_count, index="shannon") #beta diversity is a measure for similarity and overlap between samples of distributions
alpha_diversity_ec = fisher.alpha(ec_matrix_count) #alpha diversity is average diversity within community
richness_ec = specnumber(ec_matrix_count) #species richness is simple species count

timepoints = c("Day 3, morning", "Day 1, morning", "Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon",
               "Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon", "Day 1, morning", "Day 1, afternoon",
               "Day 3, morning", "Day 3, afternoon", "Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")
diversity_df_pfam = data.frame(subject, timepoints, beta_diversity_pfam, alpha_diversity_pfam, richness_pfam)
diversity_df_pfam = diversity_df_pfam[-c(1,2),]
diversity_df_ec = data.frame(subject, timepoints, beta_diversity_ec, alpha_diversity_ec, richness_ec)
diversity_df_ec = diversity_df_ec[-c(1,2),]  

beta_diversity_pfam_plot = ggplot(data = diversity_df_pfam, aes(x=timepoints, y=beta_diversity_pfam, group=subject)) +
                           geom_line(aes(color=subject), size=2) +
                           geom_point(aes(color=subject), size=4) +
                           scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                           scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                           xlab("Timepoints") + ylab("Beta diversity (Shannon Index)") +
                           labs(title="Pfam") +
                           theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                                   axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                                   legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                                   plot.title=element_text(size=20,colour="black"))
beta_diversity_pfam_plot

alpha_diversity_pfam_plot = ggplot(data = diversity_df_pfam, aes(x=timepoints, y=alpha_diversity_pfam, group=subject)) +
                            geom_line(aes(color=subject), size=2) +
                            geom_point(aes(color=subject), size=4) +
                            scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                            scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                            xlab("Timepoints") + ylab("Alpha diversity (Fisher)") +
                            labs(title="Pfam") +
                            theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                                    axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                                    legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                                    plot.title=element_text(size=20,colour="black"))
alpha_diversity_pfam_plot

richness_pfam_plot = ggplot(data = diversity_df_pfam, aes(x=timepoints, y=richness_pfam, group=subject)) +
                     geom_line(aes(color=subject), size=2) +
                     geom_point(aes(color=subject), size=4) +
                     scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                     scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                     xlab("Timepoints") + ylab("Richness") +
                     labs(title="Pfam") +
                     theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                             axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                             legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                             plot.title=element_text(size=20,colour="black"))
richness_pfam_plot

beta_diversity_ec_plot = ggplot(data = diversity_df_ec, aes(x=timepoints, y=beta_diversity_ec, group=subject)) +
                           geom_line(aes(color=subject), size=2) +
                           geom_point(aes(color=subject), size=4) +
                           scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                           scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                           xlab("Timepoints") + ylab("Beta diversity (Shannon Index)") +
                           labs(title="EC") +
                           theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                                   axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                                   legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                                   plot.title=element_text(size=20,colour="black"))
beta_diversity_ec_plot

alpha_diversity_ec_plot = ggplot(data = diversity_df_ec, aes(x=timepoints, y=alpha_diversity_ec, group=subject)) +
                            geom_line(aes(color=subject), size=2) +
                            geom_point(aes(color=subject), size=4) +
                            scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                            scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                            xlab("Timepoints") + ylab("Alpha diversity (Fisher)") +
                            labs(title="EC") +
                            theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                                    axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                                    legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                                    plot.title=element_text(size=20,colour="black"))
alpha_diversity_ec_plot

richness_ec_plot = ggplot(data = diversity_df_ec, aes(x=timepoints, y=richness_ec, group=subject)) +
                   geom_line(aes(color=subject), size=2) +
                   geom_point(aes(color=subject), size=4) +
                   scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                   scale_x_discrete(limits=c("Day 1, morning", "Day 1, afternoon", "Day 3, morning", "Day 3, afternoon")) +
                   xlab("Timepoints") + ylab("Richness") +
                   labs(title="EC") +
                   theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                           axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                           legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                           plot.title=element_text(size=20,colour="black"))
richness_ec_plot

#PCA with count data############################################################################################################################################
pfam_pca_count = prcomp(pfam_matrix_count, center=FALSE, scale=TRUE)
ec_pca_count = prcomp(ec_matrix_count, center=FALSE, scale=TRUE)

pfam_count_percentage = round(pfam_pca_count$sdev / sum(pfam_pca_count$sdev) * 100, 2)
ec_count_percentage = round(ec_pca_count$sdev / sum(ec_pca_count$sdev) * 100, 2)

pfam_count_pc = sprintf("PC%s (%s%%)", which(pfam_count_percentage==pfam_count_percentage), pfam_count_percentage)
ec_count_pc = sprintf("PC%s (%s%%)", which(ec_count_percentage==ec_count_percentage), ec_count_percentage)
  
pfam_pca_count_plot = ggplot(pfam_pca_count$x,aes(x=PC1,y=PC2,color=subject)) +
                      geom_point(size=5) +
                      xlab(pfam_count_pc[1]) + ylab(pfam_count_pc[2]) +
                      labs(title="Pfam domains") +
                      scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                      theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                              axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                              legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                              plot.title=element_text(size=20,colour="black"))
pfam_pca_count_plot

ec_pca_count_plot = ggplot(ec_pca_count$x,aes(x=PC1,y=PC2,color=subject)) +
                    geom_point(size=5) +
                    xlab(ec_count_pc[1]) + ylab(ec_count_pc[2]) +
                    labs(title="ECs") +
                    scale_color_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                    theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                            axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                            legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                            plot.title=element_text(size=20,colour="black"))
ec_pca_count_plot

#PCA with binary data###########################################################################################################################################
pfam_matrix_binary = cbind(pfam_matrix_count)
pfam_matrix_binary[pfam_matrix_binary>0] = 1

ec_matrix_binary = cbind(ec_matrix_count)
ec_matrix_binary[ec_matrix_binary>0] = 1

pfam_pca_binary = prcomp(pfam_matrix_binary, center=FALSE, scale=FALSE)
ec_pca_binary = prcomp(ec_matrix_binary, center=FALSE, scale=FALSE)

pfam_binary_percentage = round(pfam_pca_binary$sdev / sum(pfam_pca_binary$sdev) * 100, 2)
ec_binary_percentage = round(ec_pca_binary$sdev / sum(ec_pca_binary$sdev) * 100, 2)

pfam_binary_pc = sprintf("PC%s (%s%%)", which(pfam_binary_percentage==pfam_binary_percentage), pfam_binary_percentage)
ec_binary_pc = sprintf("PC%s (%s%%)", which(ec_binary_percentage==ec_binary_percentage), ec_binary_percentage)

pfam_pca_binary_plot = ggplot(pfam_pca_binary$x,aes(x=PC1,y=PC2,color=subject)) +
                      geom_point(size=5) +
                      xlab(pfam_binary_pc[1]) + ylab(pfam_binary_pc[2]) +
                      labs(title= "Pfam domains binary") +
                      scale_colour_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                      theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                              axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                              legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                              plot.title=element_text(size=20,colour="black"))
pfam_pca_binary_plot

ec_pca_binary_plot = ggplot(ec_pca_binary$x,aes(x=PC1,y=PC2,color=subject)) +
                     geom_point(size=5) +
                     xlab(ec_binary_pc[1]) + ylab(ec_binary_pc[2]) +
                     labs(title="ECs binary") +
                     scale_colour_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                     theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                             axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                             legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                             plot.title=element_text(size=20,colour="black"))
ec_pca_binary_plot

#Hierarchical clustering with count data########################################################################################################################
hc_complete_pfam = hclust(dist(pfam_matrix_count), method="complete")
hc_complete_ec = hclust(dist(ec_matrix_count), method="complete")

plot(hc_complete_pfam, main="Complete linkage, Pfam domains", xlab="", sub="", cex=1.5)
plot(hc_complete_ec, main="Complete linkage, ECs", xlab="", sub="", cex=1.5)

#Venn-diagram BETWEEN subjects##################################################################################################################################
#Obtain list of domains/enzymes in each subject
subject1_pfam = unique(pfam_count[pfam_count$Sample %in% c("NG-5450_A", "NG-5593_1A", "NG-5593_1B", "NG-5593_1C", "NG-5593_1D"),2])
subject2_pfam = unique(pfam_count[pfam_count$Sample %in% c("NG-5450_B", "NG-5593_2A", "NG-5593_2B", "NG-5593_2C", "NG-5593_2D"),2])
subject3_pfam = unique(pfam_count[pfam_count$Sample %in% c("NG-5593_3A", "NG-5593_3B", "NG-5593_3C", "NG-5593_3D"),2])
subject4_pfam = unique(pfam_count[pfam_count$Sample %in% c("NG-5593_4A", "NG-5593_4B", "NG-5593_4C", "NG-5593_4D"),2])

subject1_ec = unique(ec_count[ec_count$Sample %in% c("NG-5450_A", "NG-5593_1A", "NG-5593_1B", "NG-5593_1C", "NG-5593_1D"),2])
subject2_ec = unique(ec_count[ec_count$Sample %in% c("NG-5450_B", "NG-5593_2A", "NG-5593_2B", "NG-5593_2C", "NG-5593_2D"),2])
subject3_ec = unique(ec_count[ec_count$Sample %in% c("NG-5593_3A", "NG-5593_3B", "NG-5593_3C", "NG-5593_3D"),2])
subject4_ec = unique(ec_count[ec_count$Sample %in% c("NG-5593_4A", "NG-5593_4B", "NG-5593_4C", "NG-5593_4D"),2])

#Find domains and enzymes present or absent in all and between subjects
subject_all_pfam = Reduce(intersect, list(subject1_pfam, subject2_pfam, subject3_pfam, subject4_pfam)) #core domains
subject_1_2_pfam = intersect(subject1_pfam, subject2_pfam)
subject_1_3_pfam = intersect(subject1_pfam, subject3_pfam)
subject_1_4_pfam = intersect(subject1_pfam, subject4_pfam)
subject_2_3_pfam = intersect(subject2_pfam, subject3_pfam)
subject_2_4_pfam = intersect(subject2_pfam, subject4_pfam)
subject_3_4_pfam = intersect(subject3_pfam, subject4_pfam)
subject_1_2_3_pfam = Reduce(intersect, list(subject1_pfam, subject2_pfam, subject3_pfam))
subject_1_2_4_pfam = Reduce(intersect, list(subject1_pfam, subject2_pfam, subject4_pfam))
subject_1_3_4_pfam = Reduce(intersect, list(subject1_pfam, subject3_pfam, subject4_pfam))
subject_2_3_4_pfam = Reduce(intersect, list(subject2_pfam, subject3_pfam, subject4_pfam))

subject_all_ec = Reduce(intersect, list(subject1_ec, subject2_ec, subject3_ec, subject4_ec)) #core enzymes
subject_1_2_ec = intersect(subject1_ec, subject2_ec)
subject_1_3_ec = intersect(subject1_ec, subject3_ec)
subject_1_4_ec = intersect(subject1_ec, subject4_ec)
subject_2_3_ec = intersect(subject2_ec, subject3_ec)
subject_2_4_ec = intersect(subject2_ec, subject4_ec)
subject_3_4_ec = intersect(subject3_ec, subject4_ec)
subject_1_2_3_ec = Reduce(intersect, list(subject1_ec, subject2_ec, subject3_ec))
subject_1_2_4_ec = Reduce(intersect, list(subject1_ec, subject2_ec, subject4_ec))
subject_1_3_4_ec = Reduce(intersect, list(subject1_ec, subject3_ec, subject4_ec))
subject_2_3_4_ec = Reduce(intersect, list(subject2_ec, subject3_ec, subject4_ec))

subject1_unique_pfam = Reduce(setdiff, list(subject1_pfam, subject_1_2_pfam, subject_1_2_3_pfam, subject_1_3_pfam, subject_1_4_pfam,
                                            subject_1_2_4_pfam, subject_all_pfam, subject_1_3_4_pfam)) #domains unique for subject 1
subject2_unique_pfam = Reduce(setdiff, list(subject2_pfam, subject_1_2_pfam, subject_1_2_3_pfam, subject_2_3_pfam, subject_2_4_pfam,
                                            subject_1_2_4_pfam, subject_all_pfam, subject_2_3_4_pfam)) #domains unique for subject 2
subject3_unique_pfam = Reduce(setdiff, list(subject3_pfam, subject_1_3_pfam, subject_1_2_3_pfam, subject_2_3_pfam, subject_3_4_pfam,
                                            subject_all_pfam, subject_2_3_4_pfam)) #domains unique for subject 3
subject4_unique_pfam = Reduce(setdiff, list(subject4_pfam, subject_1_4_pfam, subject_2_4_pfam, subject_3_4_pfam, subject_1_3_4_pfam,
                                            subject_1_2_4_pfam, subject_2_3_4_pfam)) #domains unique for subject 4

subject1_unique_ec = Reduce(setdiff, list(subject1_ec, subject_1_2_ec, subject_1_2_3_ec, subject_1_3_ec, subject_1_4_ec,
                                            subject_1_2_4_ec, subject_all_ec, subject_1_3_4_ec)) #enzymes unique for subject 1
subject2_unique_ec = Reduce(setdiff, list(subject2_ec, subject_1_2_ec, subject_1_2_3_ec, subject_2_3_ec, subject_2_4_ec,
                                            subject_1_2_4_ec, subject_all_ec, subject_2_3_4_ec)) #enzymes unique for subject 2
subject3_unique_ec = Reduce(setdiff, list(subject3_ec, subject_1_3_ec, subject_1_2_3_ec, subject_2_3_ec, subject_3_4_ec,
                                            subject_all_ec, subject_2_3_4_ec)) #enzymes unique for subject 3
subject4_unique_ec = Reduce(setdiff, list(subject4_ec, subject_1_4_ec, subject_2_4_ec, subject_3_4_ec, subject_1_3_4_ec,
                                            subject_1_2_4_ec, subject_2_3_4_ec)) #enzymes unique for subject 4

#Venn-diagram of subjects
grid.newpage()
venn_pfam = draw.quad.venn(area1=length(subject1_pfam), area2=length(subject2_pfam),
               area3=length(subject3_pfam), area4=length(subject4_pfam),
               n12=length(subject_1_2_pfam), n13=length(subject_1_3_pfam),
               n14=length(subject_1_4_pfam), n23=length(subject_2_3_pfam),
               n24=length(subject_2_4_pfam), n34=length(subject_3_4_pfam),
               n1234=length(subject_all_pfam), n123=length(subject_1_2_3_pfam),
               n124=length(subject_1_2_4_pfam), n134=length(subject_1_3_4_pfam),
               n234=length(subject_2_3_4_pfam),
               category=c("Subject 1 male", "Subject 2 female", 
                          "Subject 3 male", "Subject 4 female"),
               fill=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c"),
               cex=1.5, cat.cex=1.5)
grid.text("Pfam domains", vjust=-23, gp=gpar(fontfamily="serif",cex=2))

grid.newpage()
venn_ec = draw.quad.venn(area1=length(subject1_ec), area2=length(subject2_ec),
                         area3=length(subject3_ec), area4=length(subject4_ec),
                         n12=length(subject_1_2_ec), n13=length(subject_1_3_ec),
                         n14=length(subject_1_4_ec), n23=length(subject_2_3_ec),
                         n24=length(subject_2_4_ec), n34=length(subject_3_4_ec),
                         n1234=length(subject_all_ec), n123=length(subject_1_2_3_ec),
                         n124=length(subject_1_2_4_ec), n134=length(subject_1_3_4_ec),
                         n234=length(subject_2_3_4_ec),
                         category=c("Subject 1 male", "Subject 2 female", 
                                    "Subject 3 male", "Subject 4 female"),
                         fill=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c"),
                         cex=1.5, cat.cex=1.5)
grid.text("ECs", vjust=-23, gp=gpar(fontfamily="serif",cex=2))

#Venn-diagram WITHIN subjects###################################################################################################################################
subject1_samples = c("NG-5593_1A", "NG-5593_1B", "NG-5593_1C", "NG-5593_1D")
subject2_samples = c("NG-5593_2A", "NG-5593_2B", "NG-5593_2C", "NG-5593_2D")
subject3_samples = c("NG-5593_3A", "NG-5593_3B", "NG-5593_3C", "NG-5593_3D")
subject4_samples = c("NG-5593_4A", "NG-5593_4B", "NG-5593_4C", "NG-5593_4D")

within_subjects_venn = function(subject_id, samples, pfam_ec_count) {
  t1 = unique(pfam_ec_count[pfam_ec_count$Sample == samples[1],2])
  t2 = unique(pfam_ec_count[pfam_ec_count$Sample == samples[2],2])
  t3 = unique(pfam_ec_count[pfam_ec_count$Sample == samples[3],2])
  t4 = unique(pfam_ec_count[pfam_ec_count$Sample == samples[4],2])
  t_all = Reduce(intersect, list(t1, t2, t3, t4)) #core domains
  t_1_2 = intersect(t1, t2)
  t_1_3 = intersect(t1, t3)
  t_1_4 = intersect(t1, t4)
  t_2_3 = intersect(t2, t3)
  t_2_4 = intersect(t2, t4)
  t_3_4 = intersect(t3, t4)
  t_1_2_3 = Reduce(intersect, list(t1, t2, t3))
  t_1_2_4 = Reduce(intersect, list(t1, t2, t4))
  t_1_3_4 = Reduce(intersect, list(t1, t3, t4))
  t_2_3_4= Reduce(intersect, list(t2, t3, t4))
  grid.newpage()
  venn_pfam = draw.quad.venn(area1=length(t1), area2=length(t2),
                             area3=length(t3), area4=length(t4),
                             n12=length(t_1_2), n13=length(t_1_3),
                             n14=length(t_1_4), n23=length(t_2_3),
                             n24=length(t_2_4), n34=length(t_3_4),
                             n1234=length(t_all), n123=length(t_1_2_3),
                             n124=length(t_1_2_4), n134=length(t_1_3_4),
                             n234=length(t_2_3_4),
                             category=c("Day 1, morning", "Day 1, afternoon", 
                                        "Day 3 morning", "Day 3, afternoon"),
                             fill=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c"),
                             cex=1.5, cat.cex=1.5)
  grid.text(subject_id, vjust=-23, gp=gpar(fontfamily="serif",cex=2))
}

venn_subject1_pfam = within_subjects_venn("Subject 1, male, Pfam domains", subject1_samples, pfam_count)
venn_subject2_pfam = within_subjects_venn("Subject 2, female, Pfam domains", subject2_samples, pfam_count)
venn_subject3_pfam = within_subjects_venn("Subject 3, male, Pfam domains", subject3_samples, pfam_count)
venn_subject4_pfam = within_subjects_venn("Subject 4, female, Pfam domains", subject4_samples, pfam_count)

venn_subject1_ec = within_subjects_venn("Subject 1, male, ECs", subject1_samples, ec_count)
venn_subject2_ec = within_subjects_venn("Subject 2, female, ECs", subject2_samples, ec_count)
venn_subject3_ec = within_subjects_venn("Subject 3, male, ECs", subject3_samples, ec_count)
venn_subject4_ec = within_subjects_venn("Subject 4, female, ECs", subject4_samples, ec_count)

#Rarefaction curve##############################################################################################################################################
unique_pfam = c()
for (sample in levels(pfam_count$Sample)) {
  unique_pfam[sample] = length(unique(pfam_count[pfam_count$Sample == sample,2]))
}

unique_ec = c()
for (sample in levels(ec_count$Sample)) {
  unique_ec[sample] = length(unique(ec_count[ec_count$Sample == sample,2]))
}

read_count_per_sample = read_data[,5]
sample_labels = read_data$X

rarefaction_pfam = ggplot(data=cbind(read_count_per_sample, unique_pfam), aes(x=read_count_per_sample, y=unique_pfam)) +
                   geom_point(size=5) +
                   geom_smooth(method="lm", formula=(y~log(x)), colour="#1f78b4") +
                   scale_x_continuous(expand=c(0,0), labels=scales::comma) +
                   scale_y_continuous(expand=c(0,0)) +
                   xlab("Number of reads") + ylab("Number of Pfams") + labs(title="Rarefaction curve Pfam") +
                   theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                           axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                           legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                           plot.title=element_text(size=20,colour="black"))
rarefaction_pfam

rarefaction_ec = ggplot(data=cbind(read_count_per_sample, unique_ec), aes(x=read_count_per_sample, y=unique_ec)) +
                 geom_point(size=5) +
                 geom_smooth(method="lm", formula=(y~log(x)), colour="#1f78b4") +
                 scale_x_continuous(expand=c(0,0), labels=scales::comma) +
                 scale_y_continuous(expand=c(0,0)) +
                 xlab("Number of reads") + ylab("Number of ECs") + labs(title="Rarefaction curve EC") +
                 theme_classic() + theme(axis.title.x=element_text(size=20,colour="black"), axis.text.x=element_text(size=18,colour="black"),
                                         axis.title.y=element_text(size=20,colour="black"), axis.text.y=element_text(size=18,colour="black"),
                                         legend.title=element_text(size=20,colour="black"), legend.text=element_text(size=18,colour="black"),
                                         plot.title=element_text(size=20,colour="black"))
rarefaction_ec

#Functions core domainome and enzymome##########################################################################################################################
core_domainome = pfam_matrix_count[,subject_all_pfam]
core_domainome = colMeans(core_domainome)
core_domainome = sort(core_domainome, decreasing=TRUE)

core_enzymome = ec_matrix_count[,subject_all_ec]
core_enzymome = colMeans(core_enzymome)
core_enzymome = sort(core_enzymome, decreasing=TRUE)

write.table(core_domainome, file="blastx_core_domainome.tsv", sep="\t", col.names=FALSE, quote=FALSE)
write.table(core_enzymome, file="blastx_core_enzymome.tsv", sep="\t", col.names=FALSE, quote=FALSE)

#Unique domainome and enzymome##################################################################################################################################
unique_domainome_sub1 = pfam_matrix_count[subject1_samples,subject1_unique_pfam]
unique_domainome_sub1 = colMeans(unique_domainome_sub1)
unique_domainome_sub1 = sort(unique_domainome_sub1, decreasing=TRUE)

unique_domainome_sub2 = pfam_matrix_count[subject2_samples,subject2_unique_pfam]
unique_domainome_sub2 = colMeans(unique_domainome_sub2)
unique_domainome_sub2 = sort(unique_domainome_sub2, decreasing=TRUE)

unique_domainome_sub3 = pfam_matrix_count[subject3_samples,subject3_unique_pfam]
unique_domainome_sub3 = colMeans(unique_domainome_sub3)
unique_domainome_sub3 = sort(unique_domainome_sub3, decreasing=TRUE)

unique_domainome_sub4 = pfam_matrix_count[subject4_samples,subject4_unique_pfam]
unique_domainome_sub4 = colMeans(unique_domainome_sub4)
unique_domainome_sub4 = sort(unique_domainome_sub4, decreasing=TRUE)

write.table(unique_domainome_sub1, file="blastx_unique_domainome_sub1.tsv", sep="\t", col.names=FALSE, quote=FALSE)
write.table(unique_domainome_sub2, file="blastx_unique_domainome_sub2.tsv", sep="\t", col.names=FALSE, quote=FALSE)
write.table(unique_domainome_sub3, file="blastx_unique_domainome_sub3.tsv", sep="\t", col.names=FALSE, quote=FALSE)
write.table(unique_domainome_sub4, file="blastx_unique_domainome_sub4.tsv", sep="\t", col.names=FALSE, quote=FALSE)

unique_enzymome_sub1 = ec_matrix_count[subject1_samples,subject1_unique_ec]
unique_enzymome_sub1 = colMeans(unique_enzymome_sub1)
unique_enzymome_sub1 = sort(unique_enzymome_sub1, decreasing=TRUE)

unique_enzymome_sub2 = ec_matrix_count[subject2_samples,subject2_unique_ec]
unique_enzymome_sub2 = colMeans(unique_enzymome_sub2)
unique_enzymome_sub2 = sort(unique_enzymome_sub2, decreasing=TRUE)

unique_enzymome_sub3 = ec_matrix_count[subject3_samples,subject3_unique_ec]
unique_enzymome_sub3 = colMeans(unique_enzymome_sub3)
unique_enzymome_sub3 = sort(unique_enzymome_sub3, decreasing=TRUE)

unique_enzymome_sub4 = ec_matrix_count[subject4_samples,subject4_unique_ec]
unique_enzymome_sub4 = colMeans(unique_enzymome_sub4)
unique_enzymome_sub4 = sort(unique_enzymome_sub4, decreasing=TRUE)

write.table(unique_enzymome_sub1, file="blastx_unique_enzymome_sub1.tsv", sep="\t", col.names=FALSE, quote=FALSE)
write.table(unique_enzymome_sub2, file="blastx_unique_enzymome_sub2.tsv", sep="\t", col.names=FALSE, quote=FALSE)
write.table(unique_enzymome_sub3, file="blastx_unique_enzymome_sub3.tsv", sep="\t", col.names=FALSE, quote=FALSE)
write.table(unique_enzymome_sub4, file="blastx_unique_enzymome_sub4.tsv", sep="\t", col.names=FALSE, quote=FALSE)

############
x_pfam = seq(1:50)
y_pfam = c()
for (i in x_pfam) {
  percentage = (nrow(pfam_count[pfam_count$Count < i,])/nrow(pfam_count))*100
  y_pfam = c(y_pfam, percentage)
}

df_pfam = cbind(x_pfam, y_pfam)

counts_plot_pfam = ggplot(df_pfam, aes(x=x_pfam, y=y_pfam)) +
                   geom_point() +
                   xlab("Number of reads") + ylab("Percentage of Pfam counts")
counts_plot_pfam

x_ec = seq(1:50)
y_ec = c()
for (i in x_ec) {
  percentage = (nrow(ec_count[ec_count$Count < i,])/nrow(ec_count))*100
  y_ec = c(y_ec, percentage)
}

df_ec = cbind(x_ec, y_ec)

counts_plot_ec = ggplot(df_ec, aes(x=x_ec, y=y_ec)) +
                 geom_point() +
                 xlab("Number of reads") + ylab("Percentage of EC counts")
counts_plot_ec

