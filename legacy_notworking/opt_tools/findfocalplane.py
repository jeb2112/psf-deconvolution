#!/usr/bin/env python

##########################
# findfocalplane.py
# Reads in a minc file, calculates something slice-by-slice,
# then plots the remainder using gnuplot
#
# Johnathon Walls, July 2006
##########################

from libOPT import *
import optparse
import os

usg = "Usage: %prog [options] filename1 <filename2> <filename3>"
parser = optparse.OptionParser(usage=usg)
parser.add_option("-q","--quiet",dest="verbose",default=False,action="store_false",
                  help="chatty?")
parser.add_option("-r","--range",dest="range",
                  help="What range of values does this span? (start:end:incr)")
parser.add_option("-m","--minimize_cache",dest="mincache",default=False,action="store_true",
                  help="Turn cache management on (slice by slice)")

(options, args) = parser.parse_args()
       
if(len(args)<1):
    parser.error("Incorrect arguments")

infiles = args
for fl in infiles:
    if(not os.path.exists(fl)):
        parser.error("The file %s does not exist." % fl)

for infile in infiles:
    start = 0
    end = -1
    incr = 1
    if(options.range):
        range = options.range.split(':')
        if(len(range) != 2 and len(range) != 3):
            parser.error("The range has the wrong number of elements.")
        start = float(range[0])
        end = float(range[1])
        if(start >= end):
            parser.error("The end is earlier than the start.")
        if(len(range)==3):    
            incr = float(range[2])
        if(start+incr >= end):
            parser.error("The increment does not allow for enough range.")
    
    ind = arange(start,end,incr)
    vals = calculate_from_minc(infile, method="sum", dimension="zspace", cache=options.mincache)
    
    fr = vals[len(vals)/2:1:-1]
    ls = vals[len(vals)/2:-1]
    
    rat = ls/fr
    
    if(end == -1):
        end = len(rat)
        ind = arange(start,end,incr)
    if(infile.find("ABS_psf_fdrready_fft.mnc") != -1):
        newrat = moving_average(rat,9)
        gplt.plot(ind,newrat/max(newrat))
    else:
        if(len(ind) != len(rat)):
            print "The range is not the same length as the vals, plotting anyways"
            gplt.plot(rat)
        else:
            gplt.plot(ind,rat/max(ravel(rat)))
    gplt.hold("on")
    
gplt.ytitle('ratio')
gplt.xtitle('Slice (zspace)')
gplt.title('ratio plot of %s' % (infile))