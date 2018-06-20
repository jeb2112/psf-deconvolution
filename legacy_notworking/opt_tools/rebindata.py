from scipy import *
import py_minc
import os
import libOPT

from optparse import OptionParser

usg = "Usage: %prog [options] infile outfile"

parser = OptionParser(usage = usg)
parser.add_option("-b", "--bins", dest="binning", default="2",
                  help="how many pixels to bin")
parser.add_option("-c", "--clobber", dest="clobber",action="store_true",
                  help="overwrite output file if it exists")
parser.add_option("-t","--tiff",dest="fromtiff",default=False,action="store_true",
                  help="starts from tiff files")

(options, args) = parser.parse_args()

if(len(args) < 2):
    parser.error("Incorrect arguments")
    
infile = args[0]
outfile = args[1]

if(not os.path.exists(infile)):
        parser.error("The specified input file"+ infile+ "does not exist.")

if(os.path.exists(outfile) and not options.clobber):
        parser.error("The specified output file"+ outfile+ "already exists.")

if(options.binning):
    binning = int(options.binning)

if(options.fromtiff):
    (y,x) = libOPT.get_tiff_size(infile)
else:
    (null,y,x) = libOPT.getminc_dimlengths(infile)

if(y%binning or x%binning):
    parser.error("The number of elements is not divisible by the binning of ", binning)

if(options.fromtiff):
    tf = libOPT.get_tempfile()
    libOPT.tifftominc(infile, tf)
else:
    tf = infile

libOPT.rebin_mincfile(tf, outfile, binning)
