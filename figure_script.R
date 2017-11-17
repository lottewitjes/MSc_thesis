#Script for making a histogram-like figure of cluster data from
#plantiSMASH with on the X-axis chromosomal locations, and on 
#the Y-axis the number of CD-HIT clusters and the bars colored
#by type of cluster (saccharide, terpene, alkaloid, lignan,
#putative).

#Written by Lotte Witjes.

#Set working directory and load library
dir="M:"
set.seed(50)
#install.packages("ggplot2")
library(ggplot2)
setwd(dir)

#Set the lengths of chromosomes and colors
chrs = c(30427671, (19698289+30427671), (19698289+30427671+23459830), 
         (19698289+30427671+23459830+18585056), 
         (19698289+30427671+23459830+18585056+26975502))

colors = c("#e6194b","#3cb44b","#808000","#0082c8","#f58231","#911eb4",
           "#808080","#f032e6","#008080","#800000","#000080")

#Read the data and configure
data = read.csv(file="clusters_thaliana.csv", header=TRUE)
data = data.frame(data[, c(1,2,3,4,6,8,10)])

#Make rectangles for background alternating colors for chromosomes
rects = data.frame(xstart=c(0, 30427671, (19698289+30427671), (19698289+30427671+23459830), 
                            (19698289+30427671+23459830+18585056)), 
                   xend=c(30427671, (19698289+30427671), (19698289+30427671+23459830), 
                   (19698289+30427671+23459830+18585056), 
                   (19698289+30427671+23459830+18585056+26975502)),
                   ystart=c(0,0,0,0,0),
                   yend=c(10,10,10,10,10),
                   col=c("gray70", "gray90", "gray70", "gray90", "gray70"))

#Make the plot
ggplot_data = ggplot(data, aes(x=((data[,4]-350000)), y=data[,7], fill=Type)) +
                    scale_x_continuous(name="Chromosomes", breaks=chrs, minor_breaks=NULL, labels=seq(1,5), expand=c(0,0)) +
                    scale_y_continuous(name="Number of CD-HIT clusters", breaks=seq(max(data[,7])), expand=c(0,0)) +
                    geom_histogram(stat="identity", width=(data[,6]*1000+350000)) +        #data[,6]*1044
                    geom_vline(xintercept = chrs, linetype="dotted") +
                    labs(title="BGCs in Arabidopsis thaliana TAIR10", size=15) +
                    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                    panel.background=element_blank(), axis.line=element_line(colour="black"),
                    axis.text.x=element_text(size=15, colour="black"), axis.text.y=element_text(size=15, colour="black"),
                    axis.title.x=element_text(size=15,colour="black"), axis.title.y=element_text(size=15, colour="black"),
                    plot.title=element_text(size=15, colour="black"), legend.text=element_text(size=15, colour="black"),
                    legend.title=element_text(size=15, colour="black")) +
                    scale_fill_manual(values=colors)
ggplot_data

#############################################################################################################
#Try scaling and plotting
data = read.csv(file="clusters_thaliana.csv", header=TRUE)
data = data.frame(data[, c(1,2,3,4,6,8,10)])

From.enlarged = data[,4]-50000
To.enlarged = data[,5]+50000
Size.enlarged = To.enlarged - From.enlarged
data = cbind(data, From.enlarged, To.enlarged, Size.enlarged)

ggplot_data = ggplot(data, aes(x=(data[,8]), y=data[,7], fill=Type)) +
  scale_x_continuous(name="Chromosomes", breaks=chrs, minor_breaks=NULL, labels=seq(1,5), expand=c(0,0)) +
  scale_y_continuous(name="Number of CD-HIT clusters", breaks=seq(max(data[,7])), expand=c(0,0)) +
  geom_histogram(stat="identity", width=data[,10]) +
  geom_vline(xintercept = chrs, linetype="dotted") +
  labs(title="BGCs in Arabidopsis thaliana TAIR10") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=colors)
ggplot_data

#################################################################################################################
#Try plotting per chromosome
data = read.csv(file="clusters_thaliana.csv", header=TRUE)
data = data.frame(data[, c(1,2,3,4,6,8,10)])

chrs = c(30427671, (19698289+30427671), (19698289+30427671+23459830), 
         (19698289+30427671+23459830+18585056), 
         (19698289+30427671+23459830+18585056+26975502))

colors = c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4",
           "#808080","#f032e6","#008080","#800000","#000080")

chr1 = data[which(data[,2]==1),]
chr2 = data[which(data[,2]==2),]
chr3 = data[which(data[,2]==3),]
chr4 = data[which(data[,2]==4),]
chr5 = data[which(data[,2]==5),]

#Chromosome 1
ggplot_chr1 = ggplot(chr1, aes(x=(chr1[,4])/1000000-0.035, y=chr1[,7], fill=Type)) +
              scale_x_continuous(name="Chromosomal location in Mb", expand=c(0,0), breaks=seq(0, (chrs[1]/1000000),5)) +
              scale_y_continuous(name="Number of CD-HIT clusters", breaks=seq(max(chr1[,7])), expand=c(0,0)) +
              geom_histogram(stat="identity", width=(chr1[,6]/1000+0.035)) +
              labs(title="BGCs in Arabidopsis thaliana TAIR10 chromosome 1") +
              theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.background=element_blank(), axis.line=element_line(colour="black"),
              axis.text.x=element_text(size=15, colour="black"), axis.text.y=element_text(size=15, colour="black"),
              axis.title.x=element_text(size=15,colour="black"), axis.title.y=element_text(size=15, colour="black"),
              plot.title=element_text(size=15, colour="black"), legend.text=element_text(size=15, colour="black"),
              legend.title=element_text(size=15, colour="black")) +
              scale_fill_manual(values=colors)
ggplot_chr1

#Chromosome 2
ggplot_chr2 = ggplot(chr2, aes(x=(chr2[,4]-0.035)/1000000, y=chr2[,7], fill=Type)) +
  scale_x_continuous(name="Chromosomal location in Mb", expand=c(0,0), breaks=seq(0, (chrs[2]/1000000),5)) +
  scale_y_continuous(name="Number of CD-HIT clusters", breaks=seq(max(chr2[,7])), expand=c(0,0)) +
  geom_histogram(stat="identity", width=(chr2[,6]/1000+0.035)) +
  labs(title="BGCs in Arabidopsis thaliana TAIR10 chromosome 2") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(colour="black"),
        axis.text.x=element_text(size=15, colour="black"), axis.text.y=element_text(size=15, colour="black"),
        axis.title.x=element_text(size=15,colour="black"), axis.title.y=element_text(size=15, colour="black"),
        plot.title=element_text(size=15, colour="black"), legend.text=element_text(size=15, colour="black"),
        legend.title=element_text(size=15, colour="black")) +
  scale_fill_manual(values=colors)
ggplot_chr2

#Chromosome 3
ggplot_chr3 = ggplot(chr3, aes(x=(chr3[,4]-0.035)/1000000, y=chr3[,7], fill=Type)) +
  scale_x_continuous(name="Chromosomal location in Mb", expand=c(0,0), breaks=seq(0, (chrs[3]/1000000),5)) +
  scale_y_continuous(name="Number of CD-HIT clusters", breaks=seq(max(chr3[,7])), expand=c(0,0)) +
  geom_histogram(stat="identity", width=(chr3[,6]/1000+0.035)) +
  labs(title="BGCs in Arabidopsis thaliana TAIR10 chromosome 3") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(colour="black"),
        axis.text.x=element_text(size=15, colour="black"), axis.text.y=element_text(size=15, colour="black"),
        axis.title.x=element_text(size=15,colour="black"), axis.title.y=element_text(size=15, colour="black"),
        plot.title=element_text(size=15, colour="black"), legend.text=element_text(size=15, colour="black"),
        legend.title=element_text(size=15, colour="black")) +
  scale_fill_manual(values=colors)
ggplot_chr3

#Chromosome 4
ggplot_chr4 = ggplot(chr4, aes(x=(chr4[,4]-0.035)/1000000, y=chr4[,7], fill=Type)) +
  scale_x_continuous(name="Chromosomal location in Mb", expand=c(0,0), breaks=seq(0, (chrs[4]/1000000),5)) +
  scale_y_continuous(name="Number of CD-HIT clusters", breaks=seq(max(chr4[,7])), expand=c(0,0)) +
  geom_histogram(stat="identity", width=(chr4[,6]/1000+0.035)) +
  labs(title="BGCs in Arabidopsis thaliana TAIR10 chromosome 4") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(colour="black"),
        axis.text.x=element_text(size=15, colour="black"), axis.text.y=element_text(size=15, colour="black"),
        axis.title.x=element_text(size=15,colour="black"), axis.title.y=element_text(size=15, colour="black"),
        plot.title=element_text(size=15, colour="black"), legend.text=element_text(size=15, colour="black"),
        legend.title=element_text(size=15, colour="black")) +
  scale_fill_manual(values=colors)
ggplot_chr4

#Chromosome 5
ggplot_chr5 = ggplot(chr5, aes(x=(chr5[,4]-0.035)/1000000, y=chr5[,7], fill=Type)) +
  scale_x_continuous(name="Chromosomal location in Mb", expand=c(0,0), breaks=seq(0, (chrs[5]/1000000),5)) +
  scale_y_continuous(name="Number of CD-HIT clusters", breaks=seq(max(chr5[,7])), expand=c(0,0)) +
  geom_histogram(stat="identity", width=(chr5[,6]/1000+0.035)) +
  labs(title="BGCs in Arabidopsis thaliana TAIR10 chromosome 5") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(colour="black"),
        axis.text.x=element_text(size=15, colour="black"), axis.text.y=element_text(size=15, colour="black"),
        axis.title.x=element_text(size=15,colour="black"), axis.title.y=element_text(size=15, colour="black"),
        plot.title=element_text(size=15, colour="black"), legend.text=element_text(size=15, colour="black"),
        legend.title=element_text(size=15, colour="black")) +
  scale_fill_manual(values=colors)
ggplot_chr5