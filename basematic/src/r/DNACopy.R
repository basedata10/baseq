library("DNAcopy")

CBS = function(){
  gc <- read.table(varbin.gc, header=F)
  thisRatio <- read.table(varbin.data, header=F)
  names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
  thisRatio$chrom <- chrom.numeric
  a <- thisRatio$bincount + 1
  thisRatio$ratio <- a / mean(a)
  thisRatio$gc.content <- gc$gc.content
  thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)
  
  a <- quantile(gc$bin.length, 0.985)
  thisRatioNobad <- thisRatio[which(bad[, 1] == 0), ]
  
  set.seed(25)
  
  CNA.object <- CNA(log(thisRatioNobad$lowratio, base=2), thisRatioNobad$chrom, thisRatioNobad$chrompos,
                    data.type="logratio", sampleid=sample.name)
  
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width)
  thisShort <- segment.smoothed.CNA.object[[2]]
  m <- matrix(data=0, nrow=nrow(thisRatioNobad), ncol=1)
  
  prevEnd <- 0
  
  for (i in 1:nrow(thisShort)) {
    thisStart <- prevEnd + 1
    thisEnd <- prevEnd + thisShort$num.marfk[i]
    m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
    prevEnd = thisEnd
  }
  
  thisRatioNobad$seg.mean.LOWESS <- m[, 1]
  return(list(ratioNobad=thisRatioNobad, short=thisShort))  
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