#!/usr/bin/env python

##############################
#
# correct_rolloff_scaling.py
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

usg = 'Usage: %prog [options] recon_infile rolloff_file outfile'
parser = OptionParser(usage=usg)

parser.add_option("-q","--quiet",dest="verbose",default=True,
                  action="store_false",help="chatty?")
parser.add_option("-a","--antialias",dest="antialias",default="1",
                  help="factor for view antialiasing")
parser.add_option("-s","--slice",dest="slice",default="all",
                  help="slices to correct.  default: 'all'")
parser.add_option("-c","--clobber",default=False,action="store_true",
                  help="overwrite existing files.")
parser.add_option("-f","--filter",default="0.3",
                  help="Rolloff value of the filter (default 0.3)")

(options, args) = parser.parse_args()

if(len(args)!=3):
    parser.error("The number of arguments is incorrect.")
    
(infile, rolloff_file, outfile) = args

if(not os.path.exists(infile)):
    parser.error("The infile %s does not exist." % infile)

if(not os.path.exists(rolloff_file)):
    parser.error("The rolloff file %s does not exist." % rollofffile)

if(os.path.exists(outfile) and not options.clobber):
    parser.error("That file %s already exists" % outfile)

(null, ny, nx) = libOPT.getminc_dimlengths(infile)

libOPT.correct_rolloff_scaling(infile, outfile, rolloff_file, rolloff=0.3, reconsize = (ny,nx), pjFOV=1.0, antialias=1)




