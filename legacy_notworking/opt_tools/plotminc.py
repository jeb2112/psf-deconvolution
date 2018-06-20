#!/usr/bin/env python

##########################
# plotminc.py
# Reads in a minc file, calculates something slice-by-slice,
# then plots the remainder using gnuplot
#
# Johnathon Walls, July 2006
##########################

from libOPT import *
import optparse
import os

usg = "Usage: %prog [options] filename"
parser = optparse.OptionParser(usage=usg)
parser.add_option("-m","--method",dest="method",default="sum",
                  help="What calculation to perform?  Default: sum.  Options: max, sum, std")
parser.add_option("-q","--quiet",dest="verbose",default=False,action="store_false",
                  help="chatty?")
parser.add_option("-d","--dimension",dest="dimension",default="zspace",
                  help="What dimension to do the calc along?")
parser.add_option("-r","--range",dest="range",
                  help="What range of values does this span? (start:end:incr)")
parser.add_option("-c","--cache",dest="cache",default="all",
                  help="Specify the cache size: options <all|slice>")

(options, args) = parser.parse_args()
       
if(len(args)!=1):
    parser.error("Incorrect arguments")

infile = args[0]

if(not os.path.exists(infile)):
    parser.error("The file does not exist.")

if(not options.dimension in valid_dimensions):
    parser.error("That is not a valid dimension")

if(not options.method in valid_slice_calcs):
    parser.error("That is not a valid calculation")

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
if(options.cache == "slice"):
    vals = calculate_from_minc_slicewise(infile, options.method, options.dimension)
else:
    vals = calculate_from_minc(infile, options.method, options.dimension, options.cache)
if(end == -1):
    end = len(vals)
    ind = arange(start,end,incr)

if(len(ind) != len(vals)):
    print "The range is not the same length as the vals, plotting anyways"
    plot(vals)
else:
    plot(ind,vals)
    
ylabel('%s' % options.method)
xlabel('Slice (%s)' % options.dimension)
title('%s plot of %s' % (options.method, infile))

show()
