
library(nls)

Model_S = function(x, y){
  model <- nls(y ~ a + (b/(x+c)), start=list(a=max(y), b=-5000, c=2), trace=T, 
               control = list(minFactor = 0.001, maxiter = 1000))
  fited = predict(model, list(x=new_x))
}

