#!/usr/bin/env python

##############################
#
# prepare_background.py
#
# Johnathon R. Walls, 2008
#
# TODO:
# 
#
##############################

import libOPT
from pylab import *
import os
from optparse import OptionParser
import shutil

usg = \
'''prepare_background.py converts a TIFF stack into an averaged,
single sliced MINC file
usage:  %prog [options] stackfile.tif outputfile.mnc
'''
parser = OptionParser(usage=usg)

parser.add_option("-f","--flip",dest="flip",default=False,action="store_true",
          help="rotate file by 90 degrees? Default: False")
parser.add_option("-q","--quiet",dest="verbose",action="store_false",default=True,
          help="Do not print messages to stdout.  Default: False")
parser.add_option("-c","--clobber",dest="clobber",action="store_true",default=False,
          help="Overwrite existing files.  Default: False")
parser.add_option("-l","--littleendian",dest="hack_endian",default=False,action="store_true",
          help="Hack for tiffsplit bug.  Default: False")

## Steps to prepare a background file.  BG files are assumed
## to be multi stack tiff files
## 1. Verify the file exists.
## 2. Convert the TIFF to minc
## 3. Average all slices in the minc file
## 4. Ensure that it is properly formatted for further processing

options, args = parser.parse_args()

if(len(args) != 2):
  parser.error("Incorrect number of arguments")

libOPT.verbose = options.verbose

infile = args[0]
outfile = args[1]

if(not os.path.exists(infile)):
  parser.error("The input file %s does not exist." % infile)
  
if(os.path.exists(outfile) and not options.clobber):
  parser.error("The output file %s already exists and no clobber is active." % outfile)
  
t1 = libOPT.get_tempfile()
libOPT.execute_command("python %s/tifftominc.py -c -r 0,4095 %s %s" % (libOPT.WORKDIR, infile, t1))
libOPT.execute_command("mincaverage -clobber -float -avgdim zspace %s %s" % (t1, outfile))
libOPT.execute_command("python %s/addzdim.py %s" % (libOPT.WORKDIR, outfile))

print "Done!"
