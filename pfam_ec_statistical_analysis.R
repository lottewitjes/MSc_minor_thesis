#Lotte Witjes
#lottewitjes@outlook.com
#22 of May 2018

#Load libraries
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("ggfortify")
#install.packages("VennDiagram")
#install.packages("vegan")
library(reshape2)
library(ggplot2)
library(ggfortify)
library(VennDiagram)
library(vegan)

#Set working directory
setwd("/media/lottewitjes/Lotte Witjes/MSc_minor_thesis/statistical_analysis/")

##################################reads###########################################
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
            scale_fill_discrete(name="Type of reads", breaks=c("Nonhuman.mRNA", "Human.mRNA", "rRNA", "Adaptor.and.low.quality"),
                                labels=c("Non-human mRNA", "Human mRNA", "rRNA", "Adaptor and low quality")) +
            scale_fill_brewer(palette="Paired") +
            theme_classic() + theme(axis.title.x=element_text(size=16,colour="black"), axis.text.x=element_text(size=14,colour="black"),
                                    axis.title.y=element_text(size=16,colour="black"), axis.text.y=element_text(size=14,colour="black"),
                                    legend.title=element_text(size=16,colour="black"), legend.text=element_text(size=14,colour="black"),
                                    plot.title=element_text(size=16,colour="black"))
read_plot

###############################domain/enzyme analyses and plots######################################
#Load data
pfam_count = read.table(file="pfam_count_metasapp.tsv", sep="\t", header=FALSE)
ec_count = read.table(file="ec_count_metasapp.tsv", sep="\t", header=FALSE)

#Set column names
colnames(pfam_count) = c("sample", "Pfam", "count")
colnames(ec_count) = c("sample", "EC", "count")

#Specify the labels
sample = c("m.s1.d3.mrn.rep", "f.s2.d1.mrn.rep", "m.s1.d1.mrn", "m.s1.d1.aft", "m.s1.d3.mrn", "m.s1.d3.aft",
                       "f.s2.d1.mrn", "f.s2.d1.aft", "f.s2.d3.mrn", "f.s2.d3.aft", "m.s3.d1.mrn", "m.s3.d1.aft", "m.s3.d3.mrn",
                       "m.s3.d3.aft", "f.s4.d1.mrn", "f.s4.d1.aft", "f.s4.d3.mrn", "f.s4.d3.aft")
subject = c("1 male", "2 female", "1 male", "1 male", "1 male", "1 male", "2 female", "2 female", "2 female", "2 female", "3 male", "3 male",
            "3 male", "3 male", "4 female", "4 female", "4 female", "4 female")

#Transform data into matrices
pfam_matrix_count = acast(pfam_count, sample~Pfam, value.var="count", fill=0)
ec_matrix_count = acast(ec_count, sample~EC, value.var="count", fill=0)

#Domain/enzyme alpha/beta/gamma diversity
diversity(pfam_matrix_count, index="shannon")
fisher.alpha(pfam_matrix_count)
specnumber(pfam_matrix_count)

diversity(ec_matrix_count, index="shannon") #beta diversity is a measure for similarity and overlap between samples of distributions
fisher.alpha(ec_matrix_count) #alpha diversity is average diversity within community
specnumber(ec_matrix_count) #species richness is simple species count

#Scale the count matrices
pfam_matrix_count = scale(pfam_matrix_count, center=FALSE, scale=TRUE)
ec_matrix_count = scale(ec_matrix_count, center=FALSE, scale=TRUE)

#PCA with count data
pfam_pca_count = prcomp(pfam_matrix_count, center=FALSE, scale=FALSE)
ec_pca_count = prcomp(ec_matrix_count, center=FALSE, scale=FALSE)

pfam_count_percentage = round(pfam_pca_count$sdev / sum(pfam_pca_count$sdev) * 100, 2)
ec__count_percentage = round(ec_pca_count$sdev / sum(ec_pca_count$sdev) * 100, 2)

pfam_pca_count_plot = ggplot(pfam_pca_count$x,aes(x=PC1,y=PC2,color=subject)) +
                      geom_point(size=4) +
                      xlab("PC1 (21.28%)") + ylab("PC2 (12.29%)") +
                      labs(title="Pfam domains") +
                      scale_color_brewer(palette="Paired") +
                      theme_classic() + theme(axis.title.x=element_text(size=16,colour="black"), axis.text.x=element_text(size=14,colour="black"),
                                              axis.title.y=element_text(size=16,colour="black"), axis.text.y=element_text(size=14,colour="black"),
                                              legend.title=element_text(size=16,colour="black"), legend.text=element_text(size=14,colour="black"),
                                              plot.title=element_text(size=16,colour="black"))
pfam_pca_count_plot

ec_pca_count_plot = ggplot(ec_pca_count$x,aes(x=PC1,y=PC2,color=subject)) +
                    geom_point(size=4) +
                    xlab("PC1 (24.17%)") + ylab("PC2 (11.24%)") +
                    labs(title="ECs") +
                    scale_colour_brewer(palette="Paired") +
                    theme_classic() + theme(axis.title.x=element_text(size=16,colour="black"), axis.text.x=element_text(size=14,colour="black"),
                                            axis.title.y=element_text(size=16,colour="black"), axis.text.y=element_text(size=14,colour="black"),
                                            legend.title=element_text(size=16,colour="black"), legend.text=element_text(size=14,colour="black"),
                                            plot.title=element_text(size=16,colour="black"))
ec_pca_count_plot

#PCA with binary data
pfam_matrix_binary = pfam_matrix_count
pfam_matrix_binary[pfam_matrix_binary>0] = 1

ec_matrix_binary = ec_matrix_count
ec_matrix_binary[ec_matrix_binary>0] = 1

pfam_pca_binary = prcomp(pfam_matrix_binary, center=FALSE, scale=FALSE)
ec_pca_binary = prcomp(ec_matrix_binary, center=FALSE, scale=FALSE)

pfam_binary_percentage = round(pfam_pca_binary$sdev / sum(pfam_pca_binary$sdev) * 100, 2)
ec_binary_percentage = round(ec_pca_binary$sdev / sum(ec_pca_binary$sdev) * 100, 2)

pfam_pca_binary_plot = ggplot(pfam_pca_binary$x,aes(x=PC1,y=PC2,color=subject)) +
                      geom_point(size=4) +
                      xlab("PC1 (28.97%)") + ylab("PC2 (8.70)") +
                      labs(title= "Pfam domains binary") +
                      scale_colour_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                      theme_classic() + theme(axis.title.x=element_text(size=16,colour="black"), axis.text.x=element_text(size=14,colour="black"),
                                              axis.title.y=element_text(size=16,colour="black"), axis.text.y=element_text(size=14,colour="black"),
                                              legend.title=element_text(size=16,colour="black"), legend.text=element_text(size=14,colour="black"),
                                              plot.title=element_text(size=16,colour="black"))
pfam_pca_binary_plot

ec_pca_binary_plot = ggplot(ec_pca_binary$x,aes(x=PC1,y=PC2,color=subject)) +
                     geom_point(size=4) +
                     xlab("PC1 (31.74%)") + ylab("PC2 (7.69%)") +
                     labs(title="ECs binary") +
                     scale_colour_manual(values=c("#a6cee3", "#b2df8a", "#1f78b4", "#33a02c")) +
                     theme_classic() + theme(axis.title.x=element_text(size=16,colour="black"), axis.text.x=element_text(size=14,colour="black"),
                                             axis.title.y=element_text(size=16,colour="black"), axis.text.y=element_text(size=14,colour="black"),
                                             legend.title=element_text(size=16,colour="black"), legend.text=element_text(size=14,colour="black"),
                                             plot.title=element_text(size=16,colour="black"))
ec_pca_binary_plot

#Hierarchical clustering with count data
rownames(pfam_matrix_count) = sample
hc.complete = hclust(dist(pfam_matrix_count), method="complete")
plot(hc.complete, main="Complete linkage", xlab="", sub="", cex=0.9)

#Obtain list of domains/enzymes in each subject
subject1_pfam = unique(fam_count[pfam_count$sample %in% c("NG-5450_A", "NG-5593_1A", "NG-5593_1B", "NG-5593_1C", "NG-5593_1D"),2])
subject2_pfam = unique(pfam_count[pfam_count$sample %in% c("NG-5450_B", "NG-5593_2A", "NG-5593_2B", "NG-5593_2C", "NG-5593_2D"),2])
subject3_pfam = unique(pfam_count[pfam_count$sample %in% c("NG-5593_3A", "NG-5593_3B", "NG-5593_3C", "NG-5593_3D"),2])
subject4_pfam = unique(pfam_count[pfam_count$sample %in% c("NG-5593_4A", "NG-5593_4B", "NG-5593_4C", "NG-5593_4D"),2])

subject1_ec = unique(ec_count[ec_count$sample %in% c("NG-5450_A", "NG-5593_1A", "NG-5593_1B", "NG-5593_1C", "NG-5593_1D"),2])
subject2_ec = unique(ec_count[ec_count$sample %in% c("NG-5450_B", "NG-5593_2A", "NG-5593_2B", "NG-5593_2C", "NG-5593_2D"),2])
subject3_ec = unique(ec_count[ec_count$sample %in% c("NG-5593_3A", "NG-5593_3B", "NG-5593_3C", "NG-5593_3D"),2])
subject4_ec = unique(ec_count[ec_count$sample %in% c("NG-5593_4A", "NG-5593_4B", "NG-5593_4C", "NG-5593_4D"),2])

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

#Rarefaction curve
unique_pfam = c(555,1117,528,279,1849,2020,1545,1681,2074,1804,
                3053,2150,3139,1817,2881,2791,2465,3139)
unique_ec = c(367,706,344,153,1285,1426,1141,1198,1411,1223,
              1782,1345,1778,1101,1698,1623,1395,1777)
read_count_per_sample = read_data[,5]
sample_labels = read_data$X

rarefaction_pfam = ggplot(data=cbind(read_count_per_sample, unique_pfam), aes(x=read_count_per_sample, y=unique_pfam)) +
                   geom_point(size=2) +
                   geom_smooth(method="lm", formula=(y~log(x)), colour="#1f78b4") +
                   scale_x_continuous(expand=c(0,0), labels=scales::comma) +
                   scale_y_continuous(expand=c(0,0)) +
                   xlab("Number of reads") + ylab("Number of Pfams") + labs(title="Rarefaction curve Pfam") +
                   theme_classic() + theme(axis.title.x=element_text(size=16,colour="black"), axis.text.x=element_text(size=14,colour="black"),
                                           axis.title.y=element_text(size=16,colour="black"), axis.text.y=element_text(size=14,colour="black"),
                                           legend.title=element_text(size=16,colour="black"), legend.text=element_text(size=14,colour="black"),
                                           plot.title=element_text(size=16,colour="black"))
rarefaction_pfam

rarefaction_ec = ggplot(data=cbind(read_count_per_sample, unique_ec), aes(x=read_count_per_sample, y=unique_ec)) +
                 geom_point(size=2) +
                 geom_smooth(method="lm", formula=(y~log(x)), colour="#1f78b4") +
                 scale_x_continuous(expand=c(0,0), labels=scales::comma) +
                 scale_y_continuous(expand=c(0,0)) +
                 xlab("Number of reads") + ylab("Number of ECs") + labs(title="Rarefaction curve EC") +
                 theme_classic() + theme(axis.title.x=element_text(size=16,colour="black"), axis.text.x=element_text(size=14,colour="black"),
                                         axis.title.y=element_text(size=16,colour="black"), axis.text.y=element_text(size=14,colour="black"),
                                         legend.title=element_text(size=16,colour="black"), legend.text=element_text(size=14,colour="black"),
                                         plot.title=element_text(size=16,colour="black"))
rarefaction_ec
