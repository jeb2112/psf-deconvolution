import libOPT
from optparse import OptionParser
import os
from scipy import *


usg = "Usage: %prog [options] infile outfile"
parser = OptionParser(usage=usg)

parser.add_option("-o","--outfile",dest="outfile",help="Output to a new file",default=False)
parser.add_option("-c","--clobber",dest="clobber",help="Overwrite existing file (applies to --outfile only)", )
parser.add_option("-d","--datatype",dest="datatype",help="data type to which to convert",default="short")
parser.add_option("-r","--data_range",dest="data_range",
                  help="The range of data (min,max) to keep (use single value < 1.0 for automatic histogram calculation",default="0.95")

(options,args) = parser.parse_args()

if(len(args) != 1):
    parser.error("Incorrect arguments")
    
infile = args[0]
outfile = args[1]

if(not os.path.exists(infile)):
    parser.error("The input file %s does not exist" % infile)

if(outfile):
    if(os.path.exists(outfile) and not options.clobber):
        parser.error("The specified output file %s already exists." % options.outfile)
    out = os.path.realpath(outfile)
else:
    out = os.path.realpath(infile)

a = libOPT.openmincfile(infile)
