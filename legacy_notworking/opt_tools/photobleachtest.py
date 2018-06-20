import libOPT
from scipy import *
import sys, os
infile = sys.argv[1]

if(not os.path.exists(infile)):
   print "The infile does not exist", infile

vals = libOPT.calculate_from_minc_slicewise(infile, method="sum")
p =  libOPT.exponential_fit(vals)

bins = arange(400.,typecode=float32)
ft = p[0]*exp(-p[1]*bins)
print p
gplt.plot(bins,vals,bins,ft)
