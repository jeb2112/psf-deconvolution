#!/usr/bin/env python

import libOPT
import sys, os
from optparse import OptionParser

usg = "Usage: %prog [options] infile outfile"
parser = OptionParser(usage=usg)
parser.add_option("-c","--clobber",dest="clobber",action="store_true",
                  default=False,help="Overwrite output file")
parser.add_option("-q","--quiet",dest="verbose",action="store_false",
                  default=True,help="Chatty?")
parser.add_option("-r","--range",dest="range",default="0,4095",
                  help="lower and upper valid range of the data.")
parser.add_option("-l","--littleendian",dest="hack_endian",action="store_true",
                  default=False,help="Hack to enforce little endian on tiffsplit for bugged tiff libraries")
parser.add_option("-f","--flip",dest="flip",action="store_true",default=False,
                  help="Flip the data by 90 degrees")
parser.add_option("-d","--datatype",dest="datatype",default="short",
                  help="data type to use")
parser.add_option("-s","--signed",dest="signed",default=False,action="store_true",
                  help="signed data?")
parser.add_option("-i","--invert",dest="invert",default=False,action="store_true",
                  help="invert the data?")
                  
(options, args) = parser.parse_args()

if(len(args) < 2):
    parser.error("The number of arguments is incorrect.")

infiles = args[0:-1]
outfile = args[-1]

for f in infiles:
  if(not os.path.exists(f)):
    parser.error("The input file %s does not exist." % f)

if(os.path.exists(outfile) and not options.clobber):
    parser.error("The output file %s already exists and noclobber is active." % outfile)

if(options.datatype not in libOPT.valid_minc_datatypes):
    parser.error("The data type %s is not valid." % options.datatype)
    
lower = int(options.range.split(",")[0])
upper = int(options.range.split(",")[1])
    
libOPT.verbose = options.verbose

if(len(infiles) > 1):
  if(libOPT.verbose):
    print "Converting many TIFF files"
  libOPT.manytiffs_to_minc(infiles, outfile, datatype=options.datatype, signed=options.signed, invert=options.invert, flip=options.flip, rang=(lower,upper))
  sys.exit()

infile = infiles[0] 

if(libOPT.get_tiff_numstacks(infile) > 1):
    if(libOPT.verbose):
        print "Converting a 3D file"
    libOPT.tiffstack_to_minc(infile, outfile, datatype=options.datatype, signed=options.signed, invert=options.invert, flip=options.flip, rang=(lower,upper), hack_endian=options.hack_endian)
else:
    if(libOPT.verbose):
        print "Converting a 2D file"
    libOPT.tifftominc(infile, outfile, datatype=options.datatype, signed=options.signed, invert=options.invert, flip=options.flip, rang=(lower,upper))