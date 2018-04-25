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

PEAKS = function(){
  #Could use a simple function: dist()
  amat <- matrix(data=0, nrow=1500000, ncol=1)
  counter <- 1
  for (j in 1:(len-1)) {
    for (k in (j+1):len) {
      N <- round((starts[j] - ends[j] + 1) * (starts[k] - ends[k] + 1)/1000)
      D <- abs(2^dfs$seg.mean[j] - 2^dfs$seg.mean[k])
      if (N > 0) {
        amat[(counter:(counter+N-1)), 1] <- rep.int(D, N)
        counter <- counter+N
      }
    }
  }

  a3 <- amat[(1:counter),1]
  a3.95 <- sort(a3)[round(.95*counter)]
  a3d <- density(a3[which(a3 < a3.95)], n=1000)

  if(callpeak=="callpeak"){
    cn1 <- a3d$x[which(peaks(as.vector(a3d$y), span=101))][2]
  } else {
    cn1 <- 0.5
  }
}

#Read Normalize Bin Counts
thisRatio <- read.csv(path_norm_counts, header=T)
names(thisRatio) <- c("id", "chrom", "chrompos", "abspos", "norm_count")
print(thisRatio[1:10,])

#Run the CBS Function
segments = CBS()
write.csv(segments, path_segs_out)



