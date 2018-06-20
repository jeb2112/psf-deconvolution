import libOPT
from scipy import *
from optparse import OptionParser
import os

usg = "Usage: %prog [options] filename"
parser = OptionParser(usage=usg)

parser.add_option("-s","--slice",dest="slice",default="350",
                  help="slice to use to verify stability.  (350 default)")
parser.add_option("-d","--dimension",dest="dimension",default="zspace",
                  help="The dimension along which to take the slice.")
parser.add_option("-q","--quiet",dest="verbose",default=True,action="store_false",
                  help="chatty?")

(options,args) = parser.parse_args()

libOPT.verbose = options.verbose

if(len(args) != 1):
    parser.error("Incorrect number of arguments.")

infile = args[0]

if(not os.path.exists(infile)):
    parser.error("The filename %s does not exist." % infile)

if(options.dimension not in libOPT.valid_dimensions):
    parser.error("That is not a valid dimension")

im = squeeze(libOPT.getslice(infile, int(options.slice), options.dimension))
(rows, cols) = shape(im)
    
bb = zeros((rows*2,cols),float32)

bb[0:rows,:] = im.astype(float32)
bb[rows:,:] = im.astype(float32)

imshow(bb)