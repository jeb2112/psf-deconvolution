#!/usr/bin/env python

##############################
#
# process_scan.py 
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
'''process_scan.py takes in the raw data and the background images
of OPT data and does background subtraction
usage:  %prog [options] <datafile> <backgroundfile> <output>
'''
parser = OptionParser(usage=usg)

parser.add_option("-q","--quiet",dest="verbose",action="store_false",default=True,
          help="Do not print messages to stdout.  Default: False")
parser.add_option("-m","--mode",dest="mode",default="emission",
                  help="scan mode? <transmission|emission> Default: emission.")
parser.add_option("-c","--clobber",dest="clobber",default=False,
                  action="store_true",help="overwrite files?")



options, args = parser.parse_args()

if(len(args) != 3):
  parser.error("Incorrect number of arguments")

libOPT.verbose = options.verbose

infile = args[0]
bgfile = args[1]
outfile = args[2]

if(not os.path.exists(infile)):
  parser.error("The input file %s does not exist." % infile)
  
if(not os.path.exists(bgfile)):
  parser.error("The background file %s does not exist." % bgfile)
  
if(os.path.exists(outfile) and not options.clobber):
  parser.error("The output file %s already exists and no clobber is active." % outfile)
  
if(options.mode == "emission"):
  libOPT.execute_command("python %s/slicewisemincmath.py -c -t float -d zspace -m sub -c %s %s %s" % 
                         (libOPT.WORKDIR, infile, bgfile, outfile))
elif(options.mode == "transmission"):
  libOPT.execute_command("python %s/slicewiselogtransform.py -c -t float -d zspace -c %s %s %s" % 
                         (libOPT.WORKDIR, infile, bgfile, outfile))
else:
  parser.error("The listed mode %s is not a valid mode <transmission or emission>" % options.mode)

print "Done!"
  
