#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path_bin_counts  = args[1]
path_segs = args[2]
path_out = args[3]

library("ggplot2")
library("Cairo")
library("dplyr")

plot_geonme = function(){
  
  CairoPNG(path, width=1200, height=200)
  
  p1 = ggplot(df3,aes(abspos, show))+
    geom_rect(data=subset(df_rect, c==1), inherit.aes=FALSE,
              color = "grey", fill=NA, alpha=0.5,
              aes(xmin=mm, xmax=ma, ymin=0, ymax=maxCV))+
    geom_point(size=0.2,color="dodgerblue")+
    geom_line(data=df_seg_plot,aes(x=pos,y=CN),color="red", size=0.7)+
    scale_y_continuous(limits=c(0, maxCV),expand = c(0,0))+
    scale_x_continuous(
      limits=c(1,max(df1$abspos)),
      expand = c(0,0), 
      breaks=df_label$pos,
      labels=df_label$label)+
    theme_bw()+
    ylab("Copy Number")+
    theme(
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text.x = element_text(colour="grey20",size=15),
      panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"))

      print(p1)
      dev.off()  
}

#position of genomes...
df = read.csv(path_bin_counts)


df_rect = ddply(df[,c("chrom", "abspos")], .(chrom), summarise, mm=min(abspos), ma=max(abspos));
df_rect$c = as.factor((1:dim(df_rect)[1] %% 2))

df_label = ddply(df[,c("chrom", "abspos")], .(chrom), summarise,
                 label=chrom[1],
                 pos=median(abspos)
             )

chr_counts = dim(df_rect)[1]
df_label$label[chr_counts-1] = "X"
df_label$label[chr_counts] = "Y"
df_label = df_label[seq(1,length(df_rect[,1]),2),]
print(df[1:10,])
print(df_rect)
print(df_label)


