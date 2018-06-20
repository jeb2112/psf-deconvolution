#########################
# autorecon.py
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
                  help="The offset of the rotational axis from the center of the detectors (in units of pixels."+
                  "  ('auto' for automatic determination of offset'")
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
elif (options.offset == "auto"):
    offset = calculate_offset(infile,accuracy=0.1)
    if(verbose):
        print "*** libOPT feedback ***: Proceeding with offset = ", offset
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

(views, slices, dets) = getminc_dimlengths(infile)
if(verbose):
    print "Minc dimensions are", views, slices, dets
    
if(options.slices):
    start_slice = int(options.slices.split(',')[0])
    end_slice = int(options.slices.split(',')[1])
    if(end_slice <= start_slice or end_slice > slices):
        parser.error("Slices are incorrectly formatted")
else:
    start_slice = 0
    end_slice = slices

td = make_tempdir()
if(verbose):
    print "Working in", td

## Necessary only because of the silliness of python, see the for loop for documentation
offset_str = "-o %f" % offset
pjFOV_str = "-F %f" % pjFOV
offsetview_str = "-v %d" % offsetview

if(ROI):
    ROI_str = "-r %s" % ROI
else:
    ROI_str = ""

if(reconsize):
    reconsize_str = "-r %s" % options.reconsize
else:
    reconsize_str = ""

for i in range(start_slice, end_slice):
    execute_command("mincreshape -clobber %s -start 0,%d,0 -count %d,1,%d %s %s/slice%04d.mnc" % 
                           (get_verbose_minc_option(verbose), i, views, dets, infile, td, i))
    ## Can't do this because python does not release memory, and so I run out of memory
    ## See: http://evanjones.ca/python-memory.html -- also, note that this is fixed in python 2.5
    #reconstruct("%s/slice%04d.mnc" % (td, i), "%s/recon%04d.mnc" % (td, i), offset, reconsize, pjFOV, ROI, offsetview)
    
    execute_command("python /home/jwalls/Australiawork/packages/opt_tools/individrecon.py %s %s %s %s %s %s/slice%04d.mnc %s/recon%04d.mnc" %
                    (offset_str, pjFOV_str, offsetview_str, reconsize_str, ROI_str, td, i, td, i))
    

execute_command("mincconcat %s -concat_dimension zspace %s/recon* %s" % 
                       (get_verbose_minc_option(verbose), td, outfile))

cleanup(td)
