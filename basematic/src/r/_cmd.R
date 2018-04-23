#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
source_dir = args[1]
command = args[2]

source2 = function(module){
  source(paste(source_dir, module, sep="/"), chdir = TRUE)
}

if(command == "ggplot2"){
  source2("ggplot2.R")
}

if(command == "Lowess"){
  source2("Lowess.R")
}

if(command == "CNV_cbs"){
  source2("DNACopy.R")
  print(args)
}