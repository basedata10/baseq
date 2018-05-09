#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path_norm_counts = args[1]
path_segs_out = args[2]

library("DNAcopy")

CBS = function(){

  CNA.object <- CNA(
    log(thisRatio$norm_count, base=2),
    thisRatio$chrom, 
    thisRatio$chrompos,
    data.type="logratio"
  )

  smoothed.CNA.object <- smooth.CNA(CNA.object)

  segment.smoothed.CNA.object <- segment(
      smoothed.CNA.object, 
      alpha = 0.1, 
      nperm = 1000,
      undo.splits = "sdundo", 
      undo.SD = 0.1, 
      min.width = 5
      )

  CNVsegs <- segment.smoothed.CNA.object[[2]]
  CNVsegs$CN = 2^CNVsegs$seg.mean
  return(CNVsegs)
}

#Read Normalize Bin Counts
thisRatio <- read.table(path_norm_counts, header=T)
names(thisRatio) <- c("chrom", "chrompos", "abspos", "raw", "norm_count")
print(thisRatio[1:10,])

#Run the CBS Function
print(mean(thisRatio$norm_count))
segments = CBS()
segments$CN = segments$CN*(mean(thisRatio$norm_count)/mean(segments$CN))
print(mean(segments$CN))

write.table(segments, path_segs_out, sep = "\t", quote=FALSE, row.names = FALSE)