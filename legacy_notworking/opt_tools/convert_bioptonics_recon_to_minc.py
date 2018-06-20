
#
# process_opt_data.py
# Script to convert a series of tiff files reconstructed via NRecon
# into a single 3D minc file
# Johnathon R. Walls
# July 2008

import libOPT
import sys, os
from optparse import OptionParser

usg = "Usage: %prog [options] dir outfile"
parser = OptionParser(usage=usg)
parser.add_option("-c","--clobber",dest="clobber",action="store_true",
                  default=False,help="Overwrite output file")
parser.add_option("-q","--quiet",dest="verbose",action="store_false",
                  default=True,help="Chatty?")
parser.add_option("-p","--prefix",dest="prefix",default="im_rec",
                  help="prefix of the tiff filenames")
parser.add_option("-r","--range",dest="range",default="0,65535",
                  help="lower and upper valid range of the data.")
parser.add_option("-d","--datatype",dest="datatype",default="short",
                  help="data type to use")
parser.add_option("-s","--signed",dest="signed",default=False,action="store_true",
                  help="signed data?")

(opts,args) = parser.parse_args()

if(len(args) != 2):
  parser.error("You must specify 2 arguments.")
  
indir = args[0]
outfile = args[1]

if(os.path.exists(outfile) and not opts.clobber):
  parser.error("The outfile %s already exists and clobber is not active." % outfile)
  
if(not os.path.exists(indir)):
  parser.error("The input directory %s does not exist." % indir)

lower = int(opts.range.split(",")[0])
upper = int(opts.range.split(",")[1])

logfile = "%s/%s.log" % (indir, opts.prefix)

if(not os.path.exists(logfile)):
  parser.error("The log file %s was not found." % logfile)

scan_info = libOPT.parse_bioptonics_log(logfile)

if(scan_info is False):
  parser.error("Could not parse the log file $s" % logfile)

first_slice = int(scan_info["First Section"])
last_slice = int(scan_info["Last Section"])

if(last_slice < first_slice):
  parser.error("Could not parse sections in log file.  1st: %d, last: %d" % (first_slice, last_slice))

infls = []
for i in range(first_slice, last_slice+1):
  fname = "%s/%s%04d.tif" % (indir, opts.prefix, i)
  if(not os.path.exists(fname)):
     parser.error("Could not find the file %s, aborting ..." % fname)
  infls.append(fname)

if(scan_info["Depth (bits)"] == "16"):
  stype = "short"
  ssigned = False
else:
  stype = "byte"
  ssigned = False

libOPT.manytiffs_to_minc(infls, outfile, datatype=stype, signed=ssigned, rang=(lower,upper))
libOPT.execute_command("minc_modify_header -dinsert zspace:step=%s %s" % (str(float(scan_info["Pixel Size (um)"])/1000.), outfile))
libOPT.execute_command("minc_modify_header -dinsert yspace:step=%s %s" % (str(float(scan_info["Pixel Size (um)"])/1000.), outfile))
libOPT.execute_command("minc_modify_header -dinsert xspace:step=%s %s" % (str(float(scan_info["Pixel Size (um)"])/1000.), outfile))

