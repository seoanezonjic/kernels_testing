#! /usr/bin/env Rscript


commute_time_kernel <- function(matrix){
#' Apply Commute-Time kernel opertaion to obtain a given matrix kernel
#' @param matrix to be transformed
#' @return Conmmute-Time kernel of the given matrix
#' @import MASS package(s)
#' @see Commute time kernel (active). J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
	require(MASS)
	# Obtain collapsed diagonal
	d <- rowSums(matrix)
	# Create diagonal matrix and change
	D <- diag(d) - matrix
	# Moore-Penrose pseudoinverse
	MPseudo <- ginv(D)
	# Return
	return(MPseudo)
}


laplacian_diffusion_kernel <- function(matrix, beta = 0.02){
#' Apply Difussion Kernel (Laplacian Exponential) operation to obtain a given matrix kernel
#' @param matrix to be trabsformed
#' @param beta 
#' @return Diffusion kernel of the given matrix
#' @import Matrix package(s)
#' @see Exponential Laplacian diffusion kernel(active). F Fouss 2012 | doi: 10.1016/j.neunet.2012.03.001
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
	require(Matrix)
	# Obtain collapsed diagonal
	d <- rowSums(matrix)
	# Create diagonal matrix and change
	D <- diag(d) - matrix
	# Apply exponential
	KExpo <- as.matrix(expm(-beta*D))
	# Return
	return(KExpo)
}


kernelized_adjacency_matrix <- function(matrix){
#' Apply Adjacency matrix kernelization opertaion to obtain a given matrix kernel
#' @param matrix to be transformed
#' @return kernelized adjacency matrix of the given matrix
#' @see Kernelized adjacency matrix (active). J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
	# Obtain lambda
	lambda = min(eigen(matrix)$values)
	# Kernelize
	KAM = matrix + diag(ncol(matrix)) * abs(lambda)	
	# Return
	return(KAM)
}


markov_exponential_difussion_kernel <- function(matrix, beta = 0.04){
#' Apply Markov exponential diffusion kernel opertaion to obtain a given matrix kernel
#' @param matrix to be transformed
#' @param beta 
#' @return Markov Exponential Diffusion Kernel of the given matrix
#' @import Matrix package(s)
#' @see Markov exponential diffusion kernel (active). G Zampieri 2018 | doi.org/10.1186/s12859-018-2025-5 . Taken from compute_kernel script
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
	require(Matrix)
	# Obtain collapsed diagonal
	d <- rowSums(matrix)
	#(beta/N)*(N*I - D + A)
	M = (diag(ncol(matrix))*ncol(matrix) - diag(d) + matrix) * (beta/ncol(matrix))
	MEDK <- as.matrix(expm(M))
	# Return
	return(MEDK)
}


random_forest_kernel <- function(matrix){
#' Apply Random-Forest kernel opertaion to obtain a given matrix kernel
#' @param matrix to be transformed
#' @return Random-Forest kernel of the given matrix
#' @see Random forest kernel. J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
	# Obtain collapsed diagonal
	d <- rowSums(matrix)
	# Create diagonal matrix and change
	D <- diag(d) - matrix
	# Kernel
	KRF <- solve(diag(ncol(matrix)) + D)
	# Return
	return(KRF)
}


von_neumann_diffusion_kernel <- function(matrix,penalization){
#' Apply Von-Neumann difussion kernel opertaion to obtain a given matrix kernel
#' @param matrix to be transformed
#' @param penalization impact of penalization, commonly (1, 0.5 or 0.1)
#' @return Von-Neumann difussion kernel of the given matrix
#' @see Von Neumann diffusion kernel. J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
#' @author Fernando Moreno Jabato <jabato(at)uma(dot)es>
	# Calculate alpha = impact_of_penalization (1, 0.5 or 0.1) * spectral radius of A. spectral radius of A = absolute value of max eigenvalue of A
	alpha = penalization * (max(eigen(matrix)$values) ** (-1)) 
	# Kernel
	VNDK <- solve(diag(ncol(matrix)) - matrix*alpha)
	# Return
	return(VNDK)
}

