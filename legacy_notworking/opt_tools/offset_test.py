#!/usr/bin/env python

##############################
#
# offset_test.py
#
# Johnathon R. Walls, 2006
#
# TODO:
# - brightfield
#
##############################

from scipy import *
import py_minc
import libOPT
import os
from optparse import OptionParser
import shutil

#if(__name__="__main__"):
usg = "usage: %prog [options] input_file output_file"
parser = OptionParser(usage=usg)

parser.add_option("-q","--quiet",dest="verbose",action="store_false",default=True,
				  help="Do not print messages to stdout")
parser.add_option("-c","--clobber",dest="clobber",action="store_true",default=False,
				  help="Overwrite existing files")
parser.add_option("-r","--range",metavar="RANGE",dest="range",default="-40:40",
				  help="Range of offsets to attempt")
#parser.add_option("-i","--increment",dest="increment",default="1.0",
#				  help="The increment between offset tests")
parser.add_option("-s","--slice",dest="slice", default="350",help="Slice to recon")
parser.add_option("-a","--antialias",dest="antialias",default="1",
				  help="factor for view antialiasing")
				  
(options, args) = parser.parse_args()

if(len(args) < 2):
	parser.error("No arguments")

libOPT.verbose = options.verbose

infile = args[0]
outfile = args[1]
if(len(options.range.split(":")) == 2):
	first = float(options.range.split(":")[0])
	last = float(options.range.split(":")[1])
	inc = 1.0
elif(len(options.range.split(":")) == 3):
	first = float(options.range.split(":")[0])
	last = float(options.range.split(":")[1])
	inc = float(options.range.split(":")[2])
else:
	parser.error("The offset range is incorrectly formatted.")
	
if(first >= last+inc):
	parser.error("The offset range is incorrectly formatted")

if(not os.path.exists(infile)):
	parser.error("Error: The specified input file does not exist")

if(os.path.exists(outfile) and not options.clobber):
	parser.error("Error: the specified output file already exists")

(nz,ny,nx) = libOPT.getminc_dimlengths(infile)
if(libOPT.verbose):
	print "Minc dimensions are", nz, ny, nx

if(int(options.slice) >= ny):
	parser.error("The specified slice %d is outside of the maximum number of slices" % int(options.slice))

## Everything verified, now let's begin
tf = libOPT.get_tempfile()
tdir = libOPT.make_tempdir()

libOPT.execute_command("mincreshape -clobber -float -start 0,%d,0 -count %d,1,%d %s %s" %
						  (int(options.slice), nz,nx, infile, tf))

maxdets = int(ceil(nx+2*max(abs(first),abs(last))))
j = 0
for i in arange(first,last,inc):
	libOPT.reconstruct(tf, "%s/offset_%04d.mnc" % (tdir, j), offset=i,reconsize=(maxdets,maxdets),antialias=int(options.antialias))
	j += 1
#	libOPT.execute_command("python /micehome/jwalls/workspace/tomography_recon/individrecon.py "\
#							  "-o %d -r %d,%d %s %s/offset%04d" %
#							  (i, nx+2*max(abs(first),abs(last)),nx+2*max(abs(first),abs(last)),
#							   tf, tdir, i+abs(first)))

## That's done, now let's concatenate it

libOPT.execute_command("mincconcat -2 -sequential %s -clobber -concat_dimension zspace %s/offset*.mnc %s" %
						  (libOPT.get_verbose_minc_option(options.verbose), tdir, outfile))

libOPT.cleanup(tf)
libOPT.cleanup(tdir)
	
