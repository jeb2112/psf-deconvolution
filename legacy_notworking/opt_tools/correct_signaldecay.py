
#!/usr/bin/env python

# correct_signaldecay.py
# calculates the sum of each view, fits that to an exponential decay function, inverts the 
# exponential decay to serve as correction factors, and then multiplies the original data
# by the correction factors

from libOPT import *
from pylab import *
import sys, os, optparse

usg = "Usage: %prog [options] filename"
parser = optparse.OptionParser(usage=usg)
parser.add_option("-q","--quiet",dest="verbose",default=False,action="store_false",
                  help="chatty?")
parser.add_option("-f","--file",dest="file",default=False,
                  help="dump the plotted values to a file")
parser.add_option("-c","--clobber",dest="clobber",default=False,action="store_true",
                  help="overwrite the output file")
parser.add_option("-p","--percentloss",dest="percentloss",default="25.0",
                  help="the percentage of signal loss necessary before correction is applied")
parser.add_option("-t","--testonly",dest="testonly",default=False,action="store_true",
                  help="test by fitting exponential, but do not correct")
                  
(options, args) = parser.parse_args()

if(options.testonly):
  if(len(args) != 1):
    parser.error("Incorrect args for test only")
else:
  if(len(args) != 2):
    parser.error('Incorrect arguments (%d).'% len(args))

infile = args[0]
if(not options.testonly):
  outfile = args[1]
                
if(not os.path.exists(infile)):
  parser.error( "The infile does not exist", infile)

if(not options.testonly and os.path.exists(outfile) and not options.clobber):
  parser.error( "The outfile already exists!" )

if(options.file):
  if(os.path.exists(options.file) and not options.clobber):
    parser.error( "The graphics file already exists!")
    

vals = calculate_from_minc_slicewise(infile, method="sum")
p =  exponential_fit(vals)

# Build graph 
bins = arange(len(vals),typecode=float32)
ft = exp(p[1]*bins)
sigloss = 100-ft[-1]/ft[0]*100
if(sigloss > float(options.percentloss)):
  sigstr = "correction required"
else:
  sigstr = "no correction required"
plot(bins,vals/p[0],bins,ft)
title('exp(-%1.5fx), sigloss ~ %2.3f -- %s' % (p[1], sigloss, sigstr))

# Save if asked, if not, show
if(options.file):
  savefig(options.file)
else:
  show()

if(sigloss > float(options.percentloss)):
  print "Sig loss of %2.2f is greater than threshold of %2.2f%%, proceeding with correction ..." % (sigloss, float(options.percentloss))
else:
  print "Sig loss of %2.2f is less than threshold of %2.2f%%, program exiting" % (sigloss, float(options.percentloss))
  sys.exit()


if(options.testonly):
  sys.exit()

#ft = p[0]*exp(p[1]*bins)
#print p
#gplt.plot(bins,vals,bins,ft)
correction_factors = 1./(exp(p[1]*bins))

mincmath_handler_extreme(infile, correction_factors, outfile)

print "Done!"