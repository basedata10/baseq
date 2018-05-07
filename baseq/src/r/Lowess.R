#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path_bin_counts = args[1]
path_GC = args[2]
path_out = args[3]

library("utils")
library("Cairo")
library("ggplot2")
library("plyr")


lowess_gc <- function(jtkx, jtky) {
      jtklow <- lowess(jtkx, log(jtky), f=0.05)
      jtkz <- approx(jtklow$x, jtklow$y, jtkx)
      return(exp(log(jtky) - jtkz$y))
}

df = read.table(path_bin_counts)
df_GC = read.table(path_GC, header = T)
df = cbind(df_GC, df)
df$counts_norm_length = df$counts/df[,6]
df$raw_count = df$counts_norm_length/mean(df$counts_norm_length)

#Aggregate the bin counts: 50kb-->500kb
print(df[1:10,])
df$id = as.integer(1:dim(df)[1]/10)
df_500kb = ddply(df, .(id), summarize, 
                 chr=chr[1],
                 pos=start[1],
                 abspos=median(absstart),
                 GC = mean(GC),
                 norm_count=median(counts)
             )

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
df_500kb$norm_count = lowess_gc(df_500kb$GC, df_500kb$norm_count)
print(path_out)
options(digits=3)
write.csv(df_500kb[,c(2,3,4,6)], path_out)

# #Plot CNV...
# figurefile = "./CNV_plot.png"
# CairoPNG(figurefile, width=1000, height=300)
# p2 = ggplot(df_500kb, aes(pos, GC))+
#   geom_point(size=0.3)+
#   ylim(0, 8)+
#   theme_bw()
# print(p2)
# dev.off()




