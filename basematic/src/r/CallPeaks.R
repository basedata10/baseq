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