from scipy import *
import py_minc, os
import libOPT
from optparse import OptionParser

usg = '%prog [options] file1 file2'
parser = OptionParser(usage=usg)

parser.add_option("-q","--quiet",dest="verbose",action="store_false",default=True,
                  help="chatty?")
parser.add_option("-b","--blur",dest="blurkernel",default="10",
                  help="What size of halfwidth for the Gaussian blurring kernel")
parser.add_option("-s","--slice",dest="slice",default="0",
                  help="which slice in the minc files to align")
parser.add_option("-d","--dimension",dest="dimension",default="zspace",
                  help="what dimension to take the slice along")

(opts, args) = parser.parse_args()

if(len(args) != 2):
    parser.error("Incorrect number of arguments")
    
if(not os.path.exists(args[0])):
    parser.error("The first file does not exist.")
    
if(not os.path.exists(args[1])):
    parser.error("The second file does not exist.")
    
    
if(not opts.dimension in libOPT.valid_dimensions):
    parser.error("That dimension is not valid: %s" % str(libOPT.valid_dimensions))
    
(nz1,ny1,nx1) = libOPT.getminc_dimlengths(args[0])
(nz2,ny2,nx2) = libOPT.getminc_dimlengths(args[1])

if(opts.dimension == "zspace"):
    if(nz1 < int(opts.slice)):
        parser.error("The requested slice is beyond the size of the first minc file")
   
    if(nz2 < int(opts.slice)):
        parser.error("The requested slice is beyond the size of the second minc file")

elif (opts.dimension == "yspace"):
    if(ny1 < int(opts.slice)):
        parser.error("The requested slice is beyond the size of the first minc file")
   
    if(ny2 < int(opts.slice)):
        parser.error("The requested slice is beyond the size of the second minc file")

else:
    if(nx1 < int(opts.slice)):
        parser.error("The requested slice is beyond the size of the first minc file")
   
    if(nx2 < int(opts.slice)):
        parser.error("The requested slice is beyond the size of the second minc file")

im1 = libOPT.getslice(args[0],int(opts.slice),opts.dimension)
im2 = libOPT.getslice(args[1],int(opts.slice),opts.dimension)

print libOPT.correlate_images(im1,im2,blur=int(opts.blurkernel))
