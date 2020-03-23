#! /usr/bin/env Rscript

#' @description 
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
#' @import optparse, R.matlab

##############################################################################
##                           CONFIGURE PROGRAM                              ##
##############################################################################

# Load necessary packages
suppressPackageStartupMessages(require(optparse))        # Parse script inputs
suppressPackageStartupMessages(require(R.matlab))

# Prepare input commands
option_list <- list(
  make_option(c("-f", "--file"), action="store", type="character",
              dest="file", help="Matrix File"),
  make_option(c("-t", "--text"), action="store_true", type="logical", default = FALSE,
              dest="text", help="Flag to change status to 'Plain to Matlab'"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file basename")
)

opt <- parse_args(OptionParser(option_list=option_list))


##############################################################################
##                           LOAD SOURCE DATA                               ##
##############################################################################

if(opt$text){
	# Load matrix
	source_matrix <- as.matrix(read.table(file=opt$file,sep="\t",quote=""))
}else{
	# Load matrix
	source_matrix <- readMat(opt$file)
	source_matrix <- source_matrix[['M']]
}


##############################################################################
##                               EXPORT DATA                                ##
##############################################################################

if(opt$text){
	writeMat(paste0(opt$output,".mat"),M=source_matrix)
}else{
	write.table(source_matrix, file=paste0(opt$output,".tab"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}
