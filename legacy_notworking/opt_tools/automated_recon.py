#!/usr/bin/env python

##############################
#
# automated_recon.py 
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
'''automated_recon.py takes multiple directories as arguments, and for
each of these directories, sets up a process_opt_data.py job for the sge_batch
queue.
usage:  %prog [options] <dir1> <dir2> <dir3> ...
'''
parser = OptionParser(usage=usg)

parser.add_option("-q","--quiet",dest="verbose",action="store_false",default=True,
          help="Do not print messages to stdout.  Default: False")
parser.add_option("-c","--clobber",dest="clobber",default=False,
                  action="store_true",help="overwrite files?")
parser.add_option("-f","--fdr_filter",dest="fdr_filter",default=False,
                  action="store_true",help="Perform FDR-based filtering?")
parser.add_option("-b","--background",dest="background",default=False,
                  action="store_true",help="Perform background subtraction?")
parser.add_option("-B","--bioptonics",dest="bioptonics",default=False,
                  action="store_true",help="Biotponics data?")

opts,args=parser.parse_args()


if(len(args) < 1):
  parser.error("You must specify at least one directory for automated recon.")
 
libOPT.verbose = opts.verbose
if(not opts.verbose):
  verbose_str = "-q"
else:
  verbose_str = ""
if(not opts.clobber):
  clobber_str = ""
else:
  clobber_str = "-c"

if(opts.fdr_filter):
  fdr_str = "--fdr_filter"
else:
  fdr_str = ""

if(not opts.background):
  background_str = "--bg_output none"
else:
  background_str = ""

## OPT Prototype mode
if(opts.bioptonics):
  for a in args:
    a1 = os.path.realpath(a)
    if(not os.path.exists(a1)):
      print("The directory %s does not exist, skipping ..." % a1)
      continue
    cmd = "sge_batch -J autorecon python %s/process_bioptonics_data.py %s %s %s %s" % (libOPT.WORKDIR, fdr_str, verbose_str, clobber_str, a1)
    libOPT.execute_command(cmd)
## OPT Prototype mode   
else:
  for a in args:
    a1 = os.path.realpath(a)
    if(not os.path.exists(a1)):
      print("The directory %s does not exist, skipping ..." % a1) 
      continue
    if(os.path.exists("%s/ANNOTATION_A.txt" % a1)):
      cmd = "sge_batch -J autorecon python %s/process_prototype_data.py %s %s %s %s %s" %  (libOPT.WORKDIR, background_str, fdr_str, verbose_str, clobber_str, a1)
    else:
      parser.error("Could not find ANNOTATION_A.txt")
    
    libOPT.execute_command(cmd)
  
print "Done!"
