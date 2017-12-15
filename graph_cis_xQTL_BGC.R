

#Lotte Witjes
#1 Dec 2017
#lottewitjes@outlook.com

#Set working directory and load library
dir="E:/MSc thesis/rice_indica/cis_xQTLs_BGC_leaf_seed/MSUv6.1/with 30% cutoff/"
set.seed(50)
install.packages("ggplot2")
install.packages("reshape2")
install.packages("stringi")
library(ggplot2)
library(reshape2)
library(stringi)
setwd(dir)

#Read file containing the number of local_eQTLs, non_local_eQTLs, and mQTLs per BGC. ClusterID and clustertype are also captured.
data = read.table(file="statistics.txt", header=TRUE)
x_labels = read.table(file="labels.txt", header=FALSE, sep="\t")
data$chr = factor(data$chr, levels=c(1,2,3,4,5,6,7,8,9,10,11,12),ordered=TRUE)
data = data[order(data$chr),]
x_axis_labels = unique(paste(data$cluster_ID, data$chr, data$cluster_type))
data = cbind(data, xje)
colors = c("#D5D2CA","#A59D95","#D5D2CA","#A59D95","#D5D2CA","#A59D95","#D5D2CA","#A59D95","#D5D2CA", "#A59D95", "#D5D2CA", "#005172","#AA272F","#34B233")
rects = data.frame(xstart = c(0,6.5,11.5,14.5,17.5,20.5,28.5,34.5,37.5,38.5,42.5,48.5), xend = c(6.5,11.5,14.5,17.5,20.5,28.5,34.5,37.5,38.5,42.5,48.5,49.5), col= c("#D5D2CA","#A59D95","#D5D2CA","#A59D95","#D5D2CA","#A59D95","#D5D2CA","#A59D95","#D5D2CA", "#A59D95", "#D5D2CA", "#A59D95"))

#Make a histogram with GGplot
#all
ggplot_histogram = ggplot() +
                   geom_rect(data=rects, aes(xmin=xstart, xmax=xend, ymin=-Inf, ymax=Inf, fill=rects$col), alpha=0.4) +
                   geom_bar(data=data, aes(x=paste(data$cluster_ID, data$chr,data$cluster_type), y=data$frequency, fill=data$type), stat="identity") +
                   scale_x_discrete(labels=x_labels$V1, limits=x_axis_labels, expand=c(0,0)) +
                   scale_y_continuous(expand=c(0,0)) +
                   scale_fill_manual(values=colors[-9:-1],
                                     breaks=c("local_eQTLs", "mQTLs", "non_local_eQTLs")) +
                   labs(x="Cluster type", y="Number of xQTLs", fill="Type of xQTL") +
                   theme(axis.text.x=element_text(size=18, colour="black"),
                         axis.text.y=element_text(size=18, colour="black"),
                         axis.title.x=element_text(size=18, colour="black"),
                         axis.title.y=element_text(size=20, colour="black"),
                         title=element_text(size=18, colour="black"),
                         panel.background=element_blank(),
                         legend.text=element_text(size=18, colour="black"),
                         axis.line=element_line(colour="black")) +
                   coord_flip()
                    
ggplot_histogram

