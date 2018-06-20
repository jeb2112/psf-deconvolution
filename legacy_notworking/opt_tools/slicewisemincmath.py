#!/usr/bin/env python

# slicewisemincmath.py
#
# multiplies every slice in a minc file by a single-slice minc file
#
# Johnathon Walls, 2006

from scipy import *
import libOPT
import os
from optparse import OptionParser

#if(__name__="__main__"):
usg = "usage: %prog [options] file1 file2 outfile"
parser = OptionParser(usage=usg)

parser.add_option("-q","--quiet",dest="verbose",action="store_false",default=True,
                  help="Do not print messages to stdout")

parser.add_option("-c","--clobber",dest="clobber",action="store_true",default=False,
                  help="Overwrite existing files")
parser.add_option("-d","--dimension",dest="dimension",default="zspace",
                  help="What dimension to do the calc along?")              
parser.add_option("-m","--method",dest="method",default="mult",
                  help="the mincmath operation to perform")
parser.add_option("-t","--datatype",dest="datatype",default="float",
                  help="the data type to pass along to mincmath")

(options, args) = parser.parse_args()

if(len(args) < 3):
    parser.error("No arguments")

libOPT.verbose = options.verbose

infile1 = args[0]
infile2 = args[1]
outfile = args[2]

if(not os.path.exists(infile1)):
    parser.error("The file1 does not exist.")

if(not os.path.exists(infile1)):
    parser.error("The file2 does not exist.")
    
if(os.path.exists(outfile) and not options.clobber):
    parser.error("The output file already exists.")

if(not options.dimension in libOPT.valid_dimensions):
    parser.error("That is not a valid dimension")

if(not options.method in libOPT.valid_mincmath_methods):
    parser.error("That is not a valid method.")
    
if(not options.datatype in libOPT.valid_minc_datatypes):
    parser.error("That is not a valid data type.")

(nz, ny, nx) = libOPT.getminc_dimlengths(infile1)
(nz2, ny2, nx2) = libOPT.getminc_dimlengths(infile2)

if(options.dimension == "zspace"):
    if(nz2 != 1):
        parser.error("nz of the second file must equal 1.")
    if(ny != ny2 or nx != nx2):
        parser.error("The ny/nx dimensions do not match up.")
if(options.dimension == "yspace"):
    if(ny2 != 1):
        parser.error("nz of the second file must equal 1.")
    if(ny != ny2 or nz != nz2):
        parser.error("The nz/ny dimensions do not match up.")
if(options.dimension == "xspace"):
    if(nx2 != 1):
        parser.error("nz of the second file must equal 1.")
    if(nz != nz2 or nx != nx2):
        parser.error("The nz/nx dimensions do not match up.")


libOPT.mincmath_handler(infile1, infile2, outfile, options.datatype, options.method, options.dimension)

