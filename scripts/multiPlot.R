#!/usr/bin/env Rscript

print("multiPlot")

####args
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3){
  print("USAGE: multiPlot.R element name file.vcf")
  stop()
}

element <- args[1]
name <- args[2]
vcf <- args[3]

#### imports
if (!require("pacman")) install.packages("pacman")
pacman::p_load(vcfR)
pacman::p_load("SNPfiltR")
pacman::p_load(ggplot2)

#### functions

## plots state of the variant file
multiPlot <- function(element, name, vcf){
  ##vcfR
  vcfR <- read.vcfR(vcf)
  elem <- extract.gt(vcfR, element=element, as.numeric=TRUE)
  type <- "elem"
  pdf(paste0("MP_",name,"_",type,".pdf"))
  boxplot(elem, col=2:8, las=3)
  dev.off()
  ##chrom
  chrom <- create.chromR(name="multiPlot", vcf=vcfR, verbose=TRUE)
  chrom <- proc.chromR(chrom, verbose=TRUE)
  type <- "hist"
  pdf(paste0("MP_",name,"_",type,".pdf"))
  plot(chrom)
  dev.off()
  type <- "line"
  pdf(paste0("MP_",name,"_",type,".pdf"))
  chromoqc(chrom)
  dev.off()
  ## print report
  print(chrom)
  print(vcfR)
}

#### calls
multiPlot(element, name, vcf)

quit()
