#! /usr/bin/env python

"""
@author Fernando Moreno Jabato <jabato(at)uma(dot)com
"""

#########################################################################
##                               CONFIG                                ##
#########################################################################

# Load necessary packages
from optparse import OptionParser
from numpy import linalg
from scipy.linalg import expm
import numpy as np
import scipy.io as sio
import re


# Prepare parser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type = "string",
                  help="Matrix Matlab File. Matrix to be used will be stored in M named object")
parser.add_option("-m", "--matrix", dest="matrix", type = "string", default = "M",
                  help="Matrix object name into matlab matrix file. [Default: %default]")
parser.add_option("-k", "--kernel", dest="kernel", type = "string", default = "ct",
                  help="Kernel method to be applied: Commute-Time (ct), Exponential Laplacian Kernel (el) ;[Default: %default]")
parser.add_option("-o", "--output", dest="output", type = "string",
                  help="Output file basename")

(options, args) = parser.parse_args()



#########################################################################
##                           DEFINE KERNELS                            ##
#########################################################################
'''
Apply Commute-Time kernel opertaion to obtain a given matrix kernel
	@param matrix to be transformed
	@return Conmmute-Time kernel of the given matrix
	@import MASS package(s)
	@see Commute time kernel (active). J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
	@author Fernando Moreno Jabato <jabato(at)uma(dot)es>
'''
def commute_time_kernel(matrix):
	# Obtain collapsed diagonal
	d = np.sum(matrix, axis = 1) 
	# Create diagonal matrix and change
	D = np.diag(d) - matrix 
	# Moore-Penrose pseudoinverse
	MPseudo = linalg.pinv(D)
	# Return
	return MPseudo


'''
Apply Difussion Kernel (Laplacian Exponential) operation to obtain a given matrix kernel
	@param matrix to be trabsformed
	@param beta 
	@return Diffusion kernel of the given matrix
	@import Matrix package(s)
	@see Exponential Laplacian diffusion kernel(active). F Fouss 2012 | doi: 10.1016/j.neunet.2012.03.001
	@author Fernando Moreno Jabato <jabato(at)uma(dot)es>
'''
def laplacian_diffusion_kernel(matrix, beta = 0.02):
	d = np.sum(matrix, axis = 1) 
	# Create diagonal matrix and change
	D = np.diag(d) - matrix 
	# Apply exponential
	KExpo = expm(-beta*D)
	# Return
	return KExpo

'''
Apply Adjacency matrix kernelization opertaion to obtain a given matrix kernel
	@param matrix to be transformed
	@return kernelized adjacency matrix of the given matrix
	@see Kernelized adjacency matrix (active). J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
	@author Fernando Moreno Jabato <jabato(at)uma(dot)es>
'''
def kernelized_adjacency_matrix(matrix):
	# Obtain lambda
	lmbda = min(linalg.eigvals(matrix))
	# Kernelize
	KAM = matrix + np.identity(matrix.shape[0]) * abs(lmbda)	
	# Return
	return KAM

'''
Apply Markov exponential diffusion kernel opertaion to obtain a given matrix kernel
	@param matrix to be transformed
	@param beta 
	@return Markov Exponential Diffusion Kernel of the given matrix
	@import Matrix package(s)
	@see Markov exponential diffusion kernel (active). G Zampieri 2018 | doi.org/10.1186/s12859-018-2025-5 . Taken from compute_kernel script
	@author Fernando Moreno Jabato <jabato(at)uma(dot)es>
'''
def markov_exponential_difussion_kernel(matrix, beta = 0.04):
	# Obtain collapsed diagonal
	d = np.sum(matrix, axis = 1) 
	#(beta/N)*(N*I - D + A)
	M = (np.identity(matrix.shape[0])*matrix.shape[0] - np.diag(d) + matrix) * (beta/matrix.shape[0])
	MEDK = expm(M)
	# Return
	return MEDK


'''
Apply Random-Forest kernel opertaion to obtain a given matrix kernel
	@param matrix to be transformed
	@return Random-Forest kernel of the given matrix
	@see Random forest kernel. J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
	@author Fernando Moreno Jabato <jabato(at)uma(dot)es>
'''
def random_forest_kernel(matrix):
	# Obtain collapsed diagonal
	d = np.sum(matrix, axis = 1) 
	# Create diagonal matrix and change
	D = np.diag(d) - matrix 
	# Kernel
	KRF = linalg.inv(np.identity(matrix.shape[0]) + D)
	# Return
	return KRF


'''
Apply Von-Neumann difussion kernel opertaion to obtain a given matrix kernel
	@param matrix to be transformed
	@param penalization impact of penalization, commonly (1, 0.5 or 0.1)
	@return Von-Neumann difussion kernel of the given matrix
	@see Von Neumann diffusion kernel. J.-K. Heriche 2014 | doi: 10.1091/mbc.E13-04-0221
	@author Fernando Moreno Jabato <jabato(at)uma(dot)es>
'''
def von_neumann_diffusion_kernel(matrix,penalization):
	# Calculate alpha = impact_of_penalization (1, 0.5 or 0.1) * spectral radius of A. spectral radius of A = absolute value of max eigenvalue of A
	alpha = penalization * (max(linalg.eigvals(matrix)) ** (-1)) 
	# Kernel
	VNDK = linalg.inv(np.identity(matrix.shape[0]) - matrix*alpha)
	# Return
	return VNDK


##############################################################################
##                           LOAD SOURCE DATA                               ##
##############################################################################

# Load matrix
source_matrix = sio.loadmat(options.filename)
source_matrix = source_matrix[options.matrix]

##############################################################################
##                              APPLY KERNELS                               ##
##############################################################################

kernels = {'ct' : "ct" in options.kernel,
			'el' : "el" in options.kernel,
			'ka' : "ka" in options.kernel,
			'me' : "me" in options.kernel,
			'rf' : "rf" in options.kernel,
			'vn' : "vn" in options.kernel}

if(kernels['ct']):
	# Calculate Kernel
	ct = commute_time_kernel(source_matrix)
	# Store
	sio.savemat(options.output+"_ct.mat",{'M':ct})


if(kernels['el']):
	# Calculate Kernel
	el = laplacian_diffusion_kernel(source_matrix)
	# Store
	sio.savemat(options.output+"_el.mat",{'M':el})


if(kernels['ka']):
	# Calculate Kernel
	ka = kernelized_adjacency_matrix(source_matrix)
	# Store
	sio.savemat(options.output+"_ka.mat",{'M':ka})


if(kernels['me']):
	# Calculate Kernel
	me = markov_exponential_difussion_kernel(source_matrix)
	# Store
	sio.savemat(options.output+"_me.mat",{'M':me})


if(kernels['rf']):
	# Calculate Kernel
	rf = random_forest_kernel(source_matrix)
	# Store
	sio.savemat(options.output+"_rf.mat",{'M':rf})


if(kernels['vn']):
	# Obtain penalization value
	splitted = options.kernel.split(',')
	indx = np.min(np.where(['vn' in x for x in splitted]))
	pen = float(re.sub("vn","",splitted[indx]))
	# Calculate Kernel
	vn = von_neumann_diffusion_kernel(source_matrix,pen)
	# Store
	sio.savemat(options.output+"_vn"+str(pen)+".mat",{'M':vn})
