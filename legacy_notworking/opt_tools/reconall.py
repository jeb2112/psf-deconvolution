#!/usr/bin/env python

#########################
# reconall.py
# Takes a multislice minc file projection
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
				  help="The offset of the rotational axis from the center of the detectors (in units of pixels).  Use auto for automatic offset determination.  Use twostep for the hopefully optimized auto offset calculation.")
parser.add_option("-r","--reconsize",dest="reconsize",
				  help="The reconstructed image size nx,ny")
parser.add_option("-R","--ROI",dest="ROI",
				  help="The region of interest to reconstruct.  Format:  sx,ex,sy,ey")
parser.add_option("-q","--quiet",dest="verbose",default="True",
				  action="store_false",help="Don't echo the command")
parser.add_option("-c","--clobber",dest="clobber",
				  action="store_true",help="Overwrite the destination file")
parser.add_option("-v","--offsetview",dest="offsetview",
				  help="The number of views to offset the gantry")
parser.add_option("-s","--slices",dest="slices",
				  help="The range of slices (start:end) to reconstruct")
parser.add_option("-a","--antialias",dest="antialias",default="1",
				  help="antialiasing factor along views dimension")
parser.add_option("-f","--filter",dest="filt",default="abs_hanning",
                  help="the filter to use with pjrec (e.g. abs_hanning).  see pjrec -h for valid filters")
parser.add_option("-p", "--plotfile",dest="plotfile",default=False,
                  help="the file that will contain the offset calc plot")
parser.add_option("-w","--window",dest="window",default="5.0",
                  help="The window across which to calculate various offsets.  This window is centered on the initial guess.")
parser.add_option("-A","--accuracy",dest="accuracy",default="0.2",
                  help="the accuracy with which to find the correct rotational axis")

(options, args) = parser.parse_args()

if(len(args) != 2):
	parser.error('Not enough arguments.')

infile = args[0]
outfile = args[1]

if(not os.path.exists(infile)):
	parser.error("The specified input file %s does not exist." % infile)

if(os.path.exists(outfile) and not options.clobber):
	parser.error("The specified output file %s already exists." % outfile)

if(options.plotfile and os.path.exists(options.plotfile) and not options.clobber):
  parser.error("The specified plot file %s already exists." % options.plotfile)
  
## Now that we know the files are ok, let's verify all the options

if(not options.pjFOV):
	pjFOV = 1.0
else:
	pjFOV = float(options.pjFOV)

(views, slices, dets) = getminc_dimlengths(infile)
if(options.slices):
	start_slice = int(options.slices.split(',')[0])
	end_slice = int(options.slices.split(',')[1])
	if(end_slice <= start_slice or end_slice > slices):
		parser.error("Slices are incorrectly formatted")
else:
	start_slice = 0
	end_slice = slices

if(not options.offset):
	offset = 0.
elif(options.offset == "auto"):
  offset_guess = get_coarse_offset(infile)
  offset = get_fine_offset(infile,int((end_slice-start_slice)/2+start_slice),float(options.accuracy),
               start_offset=offset_guess-float(options.window)/2.0,
               end_offset=offset_guess+float(options.window)/2.0,
               plotfile=options.plotfile)
elif(options.offset == "twostep"):
  offset_guess = get_coarse_offset(infile)
  x1 = offset_guess-float(options.window)/2.0
  x2 = offset_guess+float(options.window)/2.0
  offset_guess2 = get_fine_offset(infile,int((end_slice-start_slice)/2+start_slice),1.0,
                  start_offset=x1,
                  end_offset=x2,
                  plotfile=options.plotfile)
  offset = get_fine_offset(infile,int((end_slice-start_slice)/2+start_slice),float(options.accuracy),
               start_offset=offset_guess2-2.0, end_offset=offset_guess2+2.0,
               plotfile=options.plotfile)
else:
	offset = float(options.offset)

if(verbose):
	print "Reconstructing with offset %2.2f" % offset

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

if(verbose):
	print "Minc dimensions are", views, slices, dets
	
td = make_tempdir()
if(verbose):
	print "Working in", td

## Necessary only because of the silliness of python, see the for loop for documentation
offset_str = "-o %f" % offset
pjFOV_str = "-F %f" % pjFOV
offsetview_str = "-v %d" % offsetview
antialias_str = "-a %s" % options.antialias
filter_str = "-f %s" % options.filt
    
if(ROI):
	ROI_str = "-r %s" % ROI
else:
	ROI_str = ""

if(reconsize):
	reconsize_str = "-r %s" % options.reconsize
else:
	reconsize_str = ""

for i in range(start_slice, end_slice):
	execute_command("mincreshape -clobber %s -start 0,%d,0 -count %d,1,%d %s %s/slice%04d.mnc" % \
						   (get_verbose_minc_option(verbose), i, views, dets, infile, td, i))
	## Can't do this because python does not release memory, and so I run out of memory
	## See: http://evanjones.ca/python-memory.html -- also, note that this is fixed in python 2.5
	#reconstruct("%s/slice%04d.mnc" % (td, i), "%s/recon%04d.mnc" % (td, i), offset, reconsize, pjFOV, ROI, offsetview)
	execute_command("python /micehome/jwalls/OPT_8.04/opt_tools/individrecon.py %s %s %s %s %s %s %s %s/slice%04d.mnc %s/recon%04d.mnc" %
				    (antialias_str, offset_str, pjFOV_str, offsetview_str, reconsize_str, ROI_str, filter_str, td, i, td, i))
	

execute_command("mincconcat -2 -sequential -clobber %s -concat_dimension zspace -sequential %s/recon* %s" % 
					   (get_verbose_minc_option(verbose), td, outfile))

cleanup(td)
