#!/usr/bin/env python

##########################
# autocrop.py
# Determines the extent of data in a minc file, and crops the 
# data accordingly.
#
# Johnathon Walls, July 2006
##########################

from libOPT import *
import optparse
import os

usg = "Usage: %prog [options] filename"
parser = optparse.OptionParser(usage=usg)
parser.add_option("-q","--quiet",dest="verbose",default=False,action="store_false",
                  help="chatty?")
parser.add_option("-o","--output",dest="outfile",default=False,
                  help="file to use as output.  Otherwise does it in-place")
parser.add_option("-c","--clobber","clobber",default=False,
                  help="overwrite outputfile")
parser.add_option("-s","--cachesize",dest="cache",default="all",
                  help="Specify the cache size: options <all|slice>")
parser.add_option("-d","--dimensions",dest="dimensions",default="zspace,yspace,xspace",
                  help="comma separated list of dimensions to crop.")

(options, args) = parser.parse_args()
       
if(len(args)!=1):
    parser.error("Incorrect arguments")

infile = args[0]

if(not os.path.exists(infile)):
    parser.error("The file does not exist.")

if(os.path.exists(options.outfile) and not options.clobber):
    parser.error("The output file exists")

dims = options.dimensions.split(",")

for dim in dims:
    if(dim not in valid_dimensions):
        parser.error("%s is not a valid dimension" % dim)

sz = get_minc_dimensions(infile)
crops = {'zspace':(0,sz[0]),'yspace':(0,sz[1]),'xspace':(o,sz[2])}

for dim in dims:

    vals = calculate_from_minc(infile, 'sum', dim, options.cache)
    extent = calculate_extent(vals)
    crops[dim] = (extent[0],extent[1]-extent[0])

if(options.outfile):
    outfile = options.outfile
else:
    outfile = get_tempfile()

execute_command("mincreshape %s -start %d,%d,%d -count %d,%d,%d %s %s" % 
                (get_minc_verbose(options.verbose),
                 crops['zspace'][0],crops['yspace'][0],crops['xspace'][0],
                 crops['zspace'][1],crops['yspace'][1],crops['xspace'][1],
                 infile, outfile))

if(not options.outfile):
    execute_command("mv %s %s" % (outfile, infile))
