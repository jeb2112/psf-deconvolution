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
parser.add_option("-s","--slices",dest="slicerange",
                  help="range of slices to bin (start:end)")
parser.add_option("-q","--quiet",dest="verbose",default=False,action="store_false",
                  help="chatty?")

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

if(options.verbose):
    libOPT.verbose = options.verbose

(slices,y,x) = libOPT.getminc_dimlengths(infile)

if(options.slicerange):
    start_slice=int(options.slicerange.split(',')[0])
    end_slice=int(options.slicerange.split(',')[1])
else:
    start_slice = 0
    end_slice = slices
    
if(y%binning or x%binning):
    parser.error("The number of elements is not divisible by the binning of ", binning)

td = libOPT.make_tempdir()

for i in range(start_slice,end_slice):
    
    libOPT.execute_command("mincreshape %s -clobber -start %d,0,0 -count 1,%d,%d %s %s/slice%04d.mnc" %
                           (libOPT.get_verbose_minc_option(libOPT.verbose), i, y, x, infile, td, i))
    libOPT.execute_command("python %s/rebindata.py -b %d %s/slice%04d.mnc %s/bslice%04d.mnc" % 
                           (libOPT.WORKDIR, binning, td, i, td, i))
    #libOPT.rebin_mincfile(tf, outfile, binning)

libOPT.execute_command("mincconcat -concat_dimension zspace %s/bs* %s" % (td, outfile))
