#####################################
##### SQANTI3 report generation ######
#####################################

### Francisco J Pardo-Palacios
### Tianyuan Liu
### Last Modified: 10/04/2022 by liut33@cardiff.ac.uk

#********************** Taking arguments from python script

args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
junc.file <- args[2]
name <- args[3]
utilities.path <- args[4]
platform <- args[5]
rdata <- args[6]
busco <- args[7]
organism <- args[8]

report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
report.file <- paste(report.prefix, "Evaluation_report.html", sep="_");
bam_file <- paste(report.prefix, "corrected.bam", sep="_")


#********************** Packages (install if not found)

list_of_packages <- c("ggplot2", "scales", "knitr","rmarkdown")
req_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(req_packages)) install.packages(req_packages, repo="http://cran.rstudio.com/")

library(ggplot2)
library(scales)
library(knitr)
library(rmarkdown)
library(Rsamagents)


#********************* Run Calculation scripts

setwd(utilities.path)
source("LRGASP_calculations.R")

LRGASP_calculations_challenge3(NAME = name, out.dir = rdata,
                      class.file=class.file, junc.file=junc.file,
                      platform = platform,
                      functions.dir = utilities.path,
                      bam = bam_file,
                      organism = organism,
                      busco = busco)

#RMD = paste(utilities.path, "Evaluation_metrics.Rmd", sep = "/")

#rmarkdown::render(RMD, params = list(
#  output.directory = rdata,
#  Name = name,
#  Platform = platform ), output_file = report.file
#)