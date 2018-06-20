#!/usr/bin/env python

from libOPT import *
from pylab import *
import sys, os, optparse

usg = "Usage: %prog [options] infile outdir"
parser = optparse.OptionParser(usage=usg)
parser.add_option("-q","--quiet",dest="verbose",default=False,action="store_false",
                  help="chatty?")
parser.add_option("-c","--clobber",dest="clobber",default=False,action="store_true",
                  help="overwrite the output file")
parser.add_option("-f","--format",dest="export_format",default="pnm",
                  help="file type to export")
parser.add_option("-D","--depth",dest="export_depth",default="16",
                  help="depth (8 or 16) to use for export")
parser.add_option("-d","--dimension",dest="dimension",default="zspace",
                  help="dimension along which to write out files")
parser.add_option("-r","--range",dest="image_range",default="0,10000",
                  help="the image range to export")

(options, args) = parser.parse_args()
if(len(args) != 2):
  parser.error("Incorrect number of arguments.")

  
infile = args[0]
outdir = args[1]

infile = os.path.realpath(infile)
outdir = os.path.realpath(outdir)

if(not os.path.exists(infile)):
  parser.error("The infile %s does not exist" % infile)

if(os.path.exists(outdir) and not options.clobber):
  parser.error("The output dir %s exists and clobber is not active" % outdir)
if(os.path.exists(outdir) and not os.path.isdir(outdir)):
  parser.error("You're trying to write out to %s but that's not a directory" % outdir) 
  
if(not os.path.exists(outdir)):
  os.mkdir(outdir)

if(options.export_format not in valid_export_types):
  parser.error("The export format %s is invalid.  Try %O" % (options.export_type, valid_export_types))

if(options.export_depth != "16" and options.export_depth != "8"):
  parser.error("The export depth %s needs to be 8 or 16." % options.export_depth)
  
if(options.dimension not in valid_dimensions):
  parser.error("The dimension %s", "is not valid.  Try %O" % (options.dimension, valid_dimensions))
  
(low,high) = options.image_range.split(",")

low = float(low)
high = float(high)

if(low >= high):
  parser.error("The image range %s is incorrectly formatted." % options.image_range)
  
export_minc(infile, outdir, options.export_format, options.export_depth, options.dimension,(low,high))
  
print "Done!"