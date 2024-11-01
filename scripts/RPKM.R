#!/usr/bin/Rscript
library(optparse)
library(readr)
library(dplyr)

option_list <- list(
	make_option(c("-i","--input"), type="character", default=NULL,
              help="prefix-telescope_report.tsv",
              dest="input_filename"),
	make_option(c("-n","--number"), type="character", default=NULL,
              help="fasta number data file",
              dest="number_filename"),
	make_option(c("-o","--output_prefix"), type="character", default="out",
              help="output filename prefix [default %default]",
              dest="output_filename"),
	make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

#if(opt$v){
#  cat(paste0("PARAMETERS:\ninput (-i): ", opt$input_filename,"\n"))
#  cat(paste0("output (-o): ", opt$output_filename,"\n"))
#  cat(paste0("number (-n): ", opt$number_filename,"\n"))
#}

options(error=traceback)
parser <- OptionParser(usage = "%prog -i prefix-telescope_report.tsv -n number_filename -o output_file", option_list=option_list)
opt = parse_args(parser)

opt$output_filename = unlist(strsplit(opt$output_filename, "/"))[length(unlist(strsplit(opt$output_filename, "/")))]



num = read.table(opt$number_filename)
num2 <- num[1,]

tel <- read_tsv(opt$input_filename,
                col_names = T, skip = 1)

tel2 <- tel %>% mutate(RPKM = 
                 (1e6 * unique_count) / (transcript_length/1000 * num2)
                 ) %>%
  select(transcript, RPKM)

write_tsv(tel2, paste0(opt$output_filename, "_RPKM.tsv"))
