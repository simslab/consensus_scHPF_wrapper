#! /usr/bin/python
import numpy as np
import pandas as pd
import argparse
from numpy.random import choice
import loompy 
from scipy.sparse import csr_matrix
import multiprocessing as mp

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l','--loom',required=True,help='Path to loom.')
    parser.add_argument('-a','--attribute',required=True,help='Column attribute for downsampling.')
    parser.add_argument('-o','--output-loom',required=True,help='Path to output loom.')
    parser.add_argument('-n','--n-molecules',type=int,required=False,help='Target number of molecules for downsampling.')
    return parser

# randomly downsamples a fraction (frac) of UMI counts from a single-cell count vector
def downsample_vector(vector,frac):
    molecules = []
    N=float(sum(vector))
    M=len(vector)
    for i in range(M):
        molecules.extend([i for m in range(vector[i])])
    dsmolecules = choice(molecules,int(frac*N),replace=False)
    dsvector = [np.count_nonzero(dsmolecules==i) for i in range(M)]
    return dsvector

# randomly downsamples a fraction (frac) of UMI counts from each single-cell count vector in a count matrix
# and saves the downsampled count matrix to a file 
def downsample_matrix(matrix,frac,outfile):
    if frac<1:
        newmatrix = np.array([downsample_vector(vec,frac) for vec in matrix.T]).T
        np.savetxt(outfile,newmatrix,fmt='%d',delimiter='\t')
    else:
        np.savetxt(outfile,matrix,fmt='%d',delimiter='\t')
    return 0

parser = parse_user_input()
ui = parser.parse_args()

data=loompy.connect(ui.loom,validate=False)
avec = data.ca[ui.attribute]

matrix = data[:,:]
cts = sum(matrix)
mncts = [np.mean(cts[avec==a]) for a in set(avec)]
if not ui.n_molecules: # if user doesn't specify the number of UMI counts for downsampling, take the submatrix with minimum avg counts
    mincts = np.min(mncts)
else:
    mincts=ui.n_molecules
cols = {attr:[] for attr in data.ca.keys()}
for i,a in enumerate(set(avec)): # generate submatrix from desired column attribute (ui.attribute) and downsample
    submatrix = matrix[:,avec==a]
    frac = mincts/mncts[i]
    outfile = 'tmp_'+str(i)+'.matrix.txt'
    downsample_matrix(submatrix,frac,outfile)
    for attr in data.ca.keys(): # re-ordered column attributes for loom
        cols[attr].extend(list(data.ca[attr][avec==a]))    

matrix = None
submatrix = None

for i,a in enumerate(set(avec)): # generate merged, downsampled count matrix from downsampled submatrices in temporary files
    outfile = 'tmp_'+str(i)+'.matrix.txt'
    if i==0:
        matrix = np.loadtxt(outfile,delimiter='\t')
    else:
        matrix = np.concatenate((matrix,np.loadtxt(outfile,delimiter='\t')),axis=1)

# generate downsampled loom
rows = {attr:data.ra[attr] for attr in data.ra.keys()}
cols = {attr:np.array(cols[attr]) for attr in data.ca.keys()}
matrix = csr_matrix(matrix)
loompy.create(ui.output_loom,matrix,col_attrs=cols,row_attrs=rows)

    

