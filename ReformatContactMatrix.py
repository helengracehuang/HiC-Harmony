import sys
import numpy as np
from scipy.sparse import diags
import argparse
import os
import time
# from memory_profiler import profile

####### matrix2diag #######
# @profile
# def readContactMatrix(in_fn):
#     """
#     reads in each line in the contact matrix
#     :param in_fn the name of file to be read
#     """
#     ContactMatrix = []
#     f = open(in_fn,'r')
#     for line in f:
#         line = line.strip()
#         array = line.split()
#         ContactMatrix.append(array)
#     ContactMatrix = np.asarray(ContactMatrix)
#     return ContactMatrix

def readContactMatrix(in_fn):
    """
    reads in each line in the contact matrix
    :param in_fn the name of file to be read
    """
    ContactMatrix = np.loadtxt(in_fn)
    return ContactMatrix

def matrix2diag(contactMatrix):
    diagonals = []
    for i in range(len(contactMatrix)):
        diagonals.append(np.diag(contactMatrix, i)) # ith diagonal
        # might improve run time & conserve memory by directly outputing each line?
    return diagonals

def outputDiags(diagonals, out_fn):
    out = open(out_fn,'a+')
    for diagonal in diagonals:
        np.savetxt(out, diagonal, newline='\t', fmt='%s')
        out.write('\n')
    out.close()

####### diag2matrix #######

def readDiagonals(in_fn):
    """
    reads in each line in the contact matrix
    :param in_fn the name of file to be read
    """
    Diagonals = []
    f = open(in_fn,'r')
    for line in f:
        line = line.strip()
        array = [float(x) for x in line.split()]
        Diagonals.append(array)
    return Diagonals

def diag2matrix(diagonals):
    fullDiag = []
    fullDiag.extend(diagonals)
    fullDiag.extend(diagonals[1:])
    offsets = []
    offsets.extend(range(len(diagonals)))
    offsets.extend(reversed(range(-len(diagonals)+1, 0)))
    ContactMatrix = diags(fullDiag, offsets).toarray()
    return ContactMatrix

def outputMatrix(contactMatrix, out_fn):
    out = open(out_fn,'w+')
    for row in contactMatrix:
        s = "\t".join(map(str, row))
        out.write(s+'\n')
    out.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Reformat contact matrix for more memory efficient band-wise computation.')
    parser.add_argument('-i', '--input', required=True, dest='input_fn',
                        help='input file name for either contact matrix or diagonals')
    parser.add_argument('-o', '--output', required=True, dest='output_fn',
                        help='output file name for either contact matrix or diagonals')
    parser.add_argument('-d', '--do', required=True, dest='do',
                        help='either matrix2diag or diag2matrix')

    args = parser.parse_args()

    in_fn = args.input_fn
    out_fn = args.output_fn

    start = time.time()

    if args.do == "matrix2diag":
        contactMatrix = readContactMatrix(in_fn)
        diagonals = matrix2diag(contactMatrix)
        outputDiags(diagonals, out_fn)
    elif args.do == "diag2matrix":
        diagonals = readDiagonals(in_fn)
        contactMatrix = diag2matrix(diagonals)
        outputMatrix(contactMatrix, out_fn)
    else:
        print("Please input the correct --do flag, the only options are matrix2diag and diag2matrix")

    end = time.time()
    
    print("computation done in", end - start)
