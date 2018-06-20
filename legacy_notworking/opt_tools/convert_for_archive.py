#!/usr/bin/env python

##############################
#
# convert_for_archive.py 
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
'''convert_for_archive.py takes an OPT recon minc file, calculates its histogram,
and determines the best areas to use as min and max, and then converts the
data to 16-bit for long term storage.  It does nothing to float or byte type
data.
usage:  %prog [options] in.mnc out.mnc
'''
parser = OptionParser(usage=usg)

parser.add_option("-q","--quiet",dest="verbose",action="store_false",default=True,
          help="Do not print messages to stdout.  Default: False")
parser.add_option("-c","--clobber",dest="clobber",default=False,
                  action="store_true",help="overwrite files?")
parser.add_option("-n","--numbins",dest="numbins",default="100000",
                  help="Number of bins in histogram")
parser.add_option("-t","--threshold",dest="threshold",default="0.99998",
                  help="percent of voxels to keep.  The threshold applies on the max side")
parser.add_option("-m","--minimum",dest="minimum",default="-100",
                  help="The minimum value of the desired result.")
parser.add_option("-p","--plotfile",dest="plotfile",default=False,
                  help="The name of the file if you want to save the histogram as a graph.")
parser.add_option("-f","--histfile",dest="histfile",default=False,
                  help="The name of the file to which the histogram will be saved.")

opts,args = parser.parse_args()

if(len(args) != 2):
  parser.error("Incorrect number of arguments.")

libOPT.verbose = opts.verbose
infile = args[0]
outfile = args[1]

if(not os.path.exists(infile)): 
  parser.error("The input file %s does not exist." % infile)
if(not opts.clobber and os.path.exists(outfile)):
  parser.error("The output file %s already exists." % outfile)
if(not opts.clobber and opts.plotfile and os.path.exists(opts.plotfile)):
  parser.error("The plot file %s already exists." % opts.plotfile)
if(not opts.clobber and opts.histfile and os.path.exists(opts.histfile)):
  parser.error("The histogram file %s already exists." % opts.histfile)
  

tf = libOPT.get_tempfile()

cmd = "mincstats -min -max %s > %s" % (infile, tf)
libOPT.execute_command(cmd)
a = open(tf)
lines = a.readlines()
a.close()
libOPT.cleanup(tf)
for line in lines:
  if(line[0:3] == "Min"):
    x = line.strip().split()
    mmin = x[-1]
  elif(line[0:3] == "Max"):
    x = line.strip().split()
    mmax = x[-1]


if(opts.minimum is not False):
  min_str = "-hist_floor %s" % opts.minimum
else:
  min_str = ""
max_str = "-hist_ceil %s" % mmax

cmd = "mincstats -histogram %s -hist_bins %s %s %s %s" % (tf, opts.numbins, min_str, max_str, infile)
libOPT.execute_command(cmd)
if(opts.histfile):
  libOPT.execute_command("cp %s %s" % (tf, opts.histfile))

a = open(tf)
lines = a.readlines()
a.close()
libOPT.cleanup(tf)
bins = []
vals = []
for line in lines:
  if(line[0] == "#"):
    continue
  if(len(line) < 2):
    continue
  x = line.strip().split()
  bins.append(float(x[0]))
  vals.append(float(x[1]))
  
bins = array(bins)
vals = array(vals)
## We don't care about the first bin, especially if we've forced a vals.
minimum = vals[1:]
bins = bins[1:]
cvals = cumsum(vals)/sum(vals)

if(opts.minimum is not False):
  min_str = opts.minimum
else:
  min_str = str(bins[0])
  
maxindex = searchsorted(cvals, float(opts.threshold))
if(maxindex > len(cvals)-1):
  maxindex = len(cvals)-1
mymax = bins[maxindex]

if(opts.clobber):
  clobber_str = "-clobber"
else:
  clobber_str = ""

if(opts.plotfile):
  plot(bins,cvals,bins,vals)
  title('histogram of %s from %s to %s' % (infile, str(bins[0]), str(bins[-1])))
  xlabel('max selected: %7.5f)' % float(mymax))
  ylabel('# bins or % of total')
  savefig(opts.plotfile)

cmd = "mincmath %s -quiet -short -clamp -const2 %s %s %s %s" % (clobber_str, min_str,str(mymax),infile, outfile)
libOPT.execute_command(cmd)
                                                             
  
