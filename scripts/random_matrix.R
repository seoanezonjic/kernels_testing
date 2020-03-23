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
  make_option(c("-s", "--size"), action="store", type="numeric",
              dest="size", help="N Matrix size (NxN)"),
  make_option(c("-t", "--type"), action="store",type="character", default = "C",
              dest="type", help="Matrix type: Connectivity (C), Integer symmetric (I), Double symmetric (D) ;[Default: %default]"),
  make_option(c("-m", "--min"), action="store", type="numeric", default = 0, 
              dest="min", help="Min value to be generated. Only used in not Connectivity mode. [Default: %default]"),
  make_option(c("-M", "--max"), action="store", type="numeric", default = 100,
              dest="max", help="Max value to be generated. Only used in not Connectivity mode. [Default: %default]"),
  make_option(c("-d", "--diag"), action="store", type="numeric", default = 1,
              dest="diag", help="Diagonal value. [Default: %default]"),
  make_option(c("-r", "--readable"), action="store_true", type="logical", default = FALSE,
              dest="readable", help="Generate also a plain text matrix file"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file basename")
)

opt <- parse_args(OptionParser(option_list=option_list))

##############################################################################
##                                 RANDOMIZE                                ##
##############################################################################

# Check matrix type
mtype<- list(Conn   = grepl("C", opt$type),
             IntSim = grepl("I", opt$type),
             DouSim = grepl("D", opt$type))

# Check
if(!any(unlist(mtype))){
	stop(paste0("Matrix type not allowed (",opt$type,")"))
}

# Generate matrix
M <- matrix(rep(0,opt$size*opt$size),ncol = opt$size)

# Symmetric matrices
if(mtype$Conn | mtype$IntSim | mtype$DouSim){
	num_items <- ((opt$size*opt$size)-opt$size)/2
	# Generate data
	if(mtype$Conn){
		items <- sample(c(0,1), size = num_items, replace = T)
	}else if(mtype$IntSim){
		items <- sample(c(opt$min,opt$max), size = num_items, replace = T)
	}else if(mtype$DouSim){
		items <- runif(num_items, opt$min, opt$max)
	}

	# Store items
	M[lower.tri(M)] <- items
	M[upper.tri(M)] <- items
	diag(M) <- rep(opt$diag, opt$size)
}



##############################################################################
##                                   WRITE                                  ##
##############################################################################

# Store matrix in Matlab format
writeMat(con = opt$output,M = M)

if(opt$readable){
	write.table(M, file = paste0(opt$output,".tab"), col.names = FALSE, row.names = FALSE, sep = "\t")
}
