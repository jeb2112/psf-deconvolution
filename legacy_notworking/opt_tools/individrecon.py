#########################
# individrecon.py
# Takes a minc file projection
# plugs it through reconstruction
#########################

from scipy import *
import tempfile, os, sys, optparse, py_minc
from libOPT import *


##################################################################################
####### SCRIPT STARTS HERE #######################################################
##################################################################################


## Necessary arguments are
## inputfile
## outputfile
## OPTIONAL args are
## FOV (default 1.)
## offset (0.)
## reconsize
## ROI (default full)
## quiet (default no)
## 
usg = 'Usage: %prog [options] input_file output_file'
parser=optparse.OptionParser(usage=usg)
parser.add_option("-F","--pjFOV",dest="pjFOV",
				  help="The FOV of the projections")
parser.add_option("-o","--offset",dest="offset",
				  help="The offset of the rotational axis from the center of the detectors (in units of pixels")
parser.add_option("-r","--reconsize",dest="reconsize",
				  help="The reconstructed image size nx,ny")
parser.add_option("-R","--ROI",dest="ROI",
				  help="The region of interest to reconstruct.  Format:  sx,ex,sy,ey")
parser.add_option("-q","--quiet",dest="quiet",
				  action="store_true",help="Don't echo the command")
parser.add_option("-c","--clobber",dest="clobber",
				  action="store_true",help="Overwrite the destination file")
parser.add_option("-v","--offsetview",dest="offsetview",
				  help="The number of views to offset the gantry")
parser.add_option("-a","--antialias",dest="antialias",default=1,
				  help="The multiple of antialiasing along the view direction.  1 is no antialiasing (default)")
parser.add_option("-f","--filter",dest="filt",default="abs_hanning",
                  help="the filter to use with pjrec") 

(options, args) = parser.parse_args()

if(len(args) != 2):
	parser.error('Not enough arguments.')

infile = args[0]
outfile = args[1]

if(not os.path.exists(infile)):
	parser.error("The specified input file %s does not exist." % infile)

if(os.path.exists(outfile) and not options.clobber):
	parser.error("The specified output file %s already exists." % outfile)

## Now that we know the files are ok, let's verify all the options

if(not options.pjFOV):
	pjFOV = 1.0
else:
	pjFOV = float(options.pjFOV)

if(not options.offset):
	offset = 0.
else:
	offset = float(options.offset)

if(not options.offsetview):
	offsetview = 0
else:
	offsetview = int(options.offsetview)

if(options.reconsize): # deal with the not case later
	reconsize = options.reconsize.split(',')
	if(len(reconsize) != 2):
		parser.error("The reconsize is incorrectly formatted.")
else:
	reconsize = False

if(options.ROI): # deal with the not case later
	ROI = options.ROI.split(',')
	if(len(ROI) != 4):
		parser.error("Incorrectly formatted ROI")

	sx = float(ROI[0])
	ex = float(ROI[1])
	sy = float(ROI[2])
	ey = float(ROI[3])
	if(sx > ey or sy > ey):
		parser.error("Incorrectly ordered ROI")
else:
	ROI = False

if(not options.filt in valid_recon_filters):
    parser.error("Filter %s is not a valid filter" % options.filt)

reconstruct(infile, outfile, offset, reconsize, pjFOV, ROI, offsetview, int(options.antialias), options.filt)

