#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("utils")
library("Cairo")
library("ggplot2")
library("plyr")

path_bin_counts = args[1]
path_GC = args[2]

lowess_gc <- function(jtkx, jtky) {
      jtklow <- lowess(jtkx, log(jtky), f=0.05)
      jtkz <- approx(jtklow$x, jtklow$y, jtkx)
      return(exp(log(jtky) - jtkz$y))
}

df = read.table(path_bin_counts)
df_GC = read.table(path_GC, header = F)
df = cbind(df_GC, df)
df$counts_norm_length = df$counts/df[,6]
df$norm_count = df$counts_norm_length/mean(df$counts_norm_length)

#Aggregate the bin counts: 50kb-->500kb
df$id = as.integer(1:dim(df)[1]/20)

#Correct by GC content
df$norm_GC = lowess_gc(df$V7, df$norm_count)*2
df_500kb = ddply(df, .(id), summarize, pos=median(V3), GC = mean(V7), norm_count=median(norm_count))

#GC Figures
figurefile = "./GC_counts.png"
CairoPNG(figurefile, width=700, height=400)
p2 = ggplot(df_500kb, aes(GC, norm_count))+
  geom_point(size=0.2)+
  ylim(0, 4)+
  xlim(0.3, 0.7)+
  theme_bw()
print(p2)
dev.off()

#Correct by GC content
df_500kb$norm_GC = lowess_gc(df_500kb$GC, df_500kb$norm_count)*2
print(df_500kb[1:100,])

#Plot CNV...
figurefile = "./CNV_plot.png"
CairoPNG(figurefile, width=1000, height=300)
p2 = ggplot(df_500kb, aes(pos, norm_GC))+
  geom_point(size=0.3)+
  ylim(0, 8)+
  theme_bw()
print(p2)
dev.off()

