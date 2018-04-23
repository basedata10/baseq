#This is code to build method for generating ...
#Code developed by Basematics
library("ggplot2")
library("Cairo")

plotFigure = function(width, height, type, path){
  if(type == "svg"){
    CairoSVG(path, width=width, height=height)
  } else {
    CairoPNG(path, width=width, height=height)
  }
  df = data.frame(x=1:10, y=(1:10)**2)
  p1 = ggplot(df, aes(x,y))+
    geom_point()+
    theme_bw()
  print(p1)
  dev.off()
}

print("Loading ggplot2 Logics...")