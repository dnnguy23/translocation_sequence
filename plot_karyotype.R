#!/usr/bin/env Rscript

options(warn=-1)

suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="Input bedGraph file"),
  make_option(c("-o", "--output"), action="store", default='bed_count/bed_cov_plot.pdf', 
              type='character', help="Output plot file"),
  make_option(c("-c", "--chromosome"), action="store", default=NA, type='character',
              help="Chromosome name"),
  make_option(c("-s", "--start"), action="store", default='0', type='character',
              help="Start coordinate"),
  make_option(c("-e", "--end"), action="store", default='0', type="character",
              help="End coordinate")
)

opt = parse_args(OptionParser(option_list=option_list))

input_file <- toString(opt$i)
output_file <- toString(opt$o)
chr_name <- opt$c
start <- as.integer(as.numeric(opt$s))
end <- as.integer(as.numeric(opt$e))

data <- read.csv(input_file, sep="\t", header=FALSE)
colnames(data) <- c('Chrom', 'Start', 'End', 'Coverage')

suppressPackageStartupMessages(library(karyoploteR))
custom.genome <- toGRanges(data.frame(chr=data$Chrom, start=data$Start, end=data$End, 
                                      name=data$Coverage))

pdf(file = output_file)
if (is.na(chr_name)) {
  kp <- plotKaryotype()
} else {
  zoom.region <- toGRanges(data.frame(chr_name, start, end))
  kp <- plotKaryotype(zoom=zoom.region)
  kpAddBaseNumbers(kp)
  kpAddCytobandLabels(kp)
}

kpBars(kp, data=custom.genome, y1=custom.genome$name, ymax=max(custom.genome$name)/4, 
       col="#003399", border="#003399")
dev.off()
