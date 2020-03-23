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
              dest="file", help="Matrix Matlab File. Matrix to be used will be stored in M named object"),
  make_option(c("-m", "--matrix"), action="store", type="character", default = "M",
              dest="matrix", help="Matrix object name into matlab matrix file. [Default: %default]"),
  make_option(c("-k", "--kernel"), action="store",type="character", default = "ct",
              dest="kernel", help="Kernel method to be applied:\n\t\t> Commute-Time (ct)\n\t\t> Exponential Laplacian Kernel (el)\n\t\t> Kernelized adjacency matrix (ka)\n\t\t> Von-Neumann Difussion Kernel with penalization * (vn*)\n\t\tNote 1: use commas (',') to separate kernel ids\n\t\tNote 2: substitute * by numbers\n\t\tNote 3: Default kernel = %default"),
  make_option(c("-o","--output"), action="store",type="character",
              dest="output", help="Output file basename")
)

opt <- parse_args(OptionParser(option_list=option_list))

full.fpath <- tryCatch(normalizePath(parent.frame(2)$ofile),  # works when using source
                       error=function(e) # works when using R CMD
                         normalizePath(unlist(strsplit(commandArgs()[grep('^--file=', commandArgs())], '='))[2]))
bname <- dirname(full.fpath)

# Import functions
if(is.null(bname)){
  source("kernels_functions.R") 
}else{
  source(paste(bname,"kernels_functions.R",sep=.Platform$file.sep))
}
to_remove <- c("to_remove","bname","source_files","script.name","script.path","initial.options")
rm(list = ls()[which(ls() %in% to_remove)])


##############################################################################
##                           LOAD SOURCE DATA                               ##
##############################################################################

# Load matrix
source_matrix <- readMat(opt$file)
source_matrix <- source_matrix[[opt$matrix]]


##############################################################################
##                              APPLY KERNELS                               ##
##############################################################################

# Check target kernels
kernels <- list(ct = grepl("ct", opt$kernel),  # Commute Time Kernel
                el = grepl("el", opt$kernel),  # Difussion Kernel (exponential lapplacian)
                ka = grepl("ka", opt$kernel),  # Kernelized Adjacency matrix
                me = grepl("me", opt$kernel),  # Markov exponential diffusion kernel
                rf = grepl("rf", opt$kernel),  # Random Forest Kernel
                vn = grepl("vn", opt$kernel)   # Von-Neumann Diffusion kernel
                )

# Check
if(!any(unlist(kernels))){
  message(paste("Kernel given is not supported :",opt$kernel))
}

# Apply kernels
if(kernels$ct){
	require(MASS)
	# Calculate Kernel
	ct <- commute_time_kernel(source_matrix)
  # Store
  writeMat(con = paste0(opt$output,"_ct.mat"),M=ct)
}
if(kernels$el){
	require(Matrix)
	# Calculate Kernel
	el <- laplacian_diffusion_kernel(source_matrix)
  # Store
  writeMat(con = paste0(opt$output,"_el.mat"),M=el)
}
if(kernels$ka){
  # Calculate Kernel
  ka <- kernelized_adjacency_matrix(source_matrix)
  # Store
  writeMat(con = paste0(opt$output,"_ka.mat"),M=ka)
}
if(kernels$me){
  # Calculate Kernel
  me <- markov_exponential_difussion_kernel(source_matrix)
  # Store
  writeMat(con = paste0(opt$output,"_me.mat"),M=me)
}
if(kernels$rf){
  # Calculate Kernel
  rf <- random_forest_kernel(source_matrix)
  # Store
  writeMat(con = paste0(opt$output,"_rf.mat"),M=rf)
}
if(kernels$vn){
  # Obtain penalization impact
  splitted <- unlist(strsplit(opt$kernel,","))
  indx <- which(grepl("vn",splitted))
  if(length(indx)<=0){
    stop("Error: you must give a number when call VN kernel")
  }
  pen <- as.numeric(gsub("vn","",splitted[indx[1]]))
  # Calculate Kernel
  vn <- von_neumann_diffusion_kernel(source_matrix,pen)
  # Store
  writeMat(con = paste0(opt$output,"_vn",gsub("vn","",splitted[indx[1]]),".mat"),M=vn)
}
