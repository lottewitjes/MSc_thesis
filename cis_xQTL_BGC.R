

#Lotte Witjes
#1 Dec 2017
#lottewitjes@outlook.com

#Set working directory and load library
dir="M:/rice_indica/cis_xQTLs_BGC_leaf_seed/"
set.seed(50)
install.packages("ggplot2")
install.packages("reshape2")
library(ggplot2)
library(reshape2)
setwd(dir)

#Read file containing the number of local_eQTLs, non_local_eQTLs, and mQTLs per BGC. ClusterID and clustertype are also captured.
data = read.table(file="statistics.txt", header=TRUE)
data$chr = factor(data$chr, levels=c(1,2,3,4,5,6,7,8,9,10,11,12),ordered=TRUE)
x_axis_labels = unique(paste(data$chr, data$cluster_type))
colors = c("#005172", "#6AADE4", "#34B233")
rects = data.frame(xstart=c(0), 
                   xend=c(4.5),
                   col=c("gray70"))

#Make a histogram with GGplot
ggplot_histogram = ggplot() +
                   geom_bar(data=data, aes(x=paste(data$chr,data$cluster_type), y=data$frequency, fill=data$type), stat="identity") +
                   scale_x_discrete(limits=x_axis_labels, expand=c(0,0)) +
                   scale_y_continuous(expand=c(0,0)) +
                   scale_fill_manual(values=colors ,
                                     breaks=c("local_eQTLs", "mQTLs", "non_local_eQTLs")) +
                   labs(title="BGCs with overlap in Oryza sativa Indica", x="Chromosome and cluster type", y="Number of xQTLs", fill="Type of xQTL") +
                   theme(axis.text.x=element_text(size=15, angle=90, hjust=1, vjust=0.5, colour="black"),
                         axis.text.y=element_text(size=15, colour="black"),
                         axis.title.x=element_text(size=15, colour="black"),
                         axis.title.y=element_text(size=15, colour="black"),
                         title=element_text(size=15, colour="black"),
                         panel.background=element_blank(),
                         legend.text=element_text(size=15, colour="black"),
                         axis.line=element_line(colour="black"))
                    
ggplot_histogram
