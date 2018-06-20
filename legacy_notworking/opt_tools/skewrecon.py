#########################
# reconall.py
# Takes a multislice minc file projection
# plugs it through reconstruction
#########################

from scipy import *
from scipy.optimize import leastsq
import tempfile, os, sys, optparse, py_minc
from libOPT import *

def straightline(parameters, y, x):
    slope, intercept= parameters
    err = y - (float(slope)*x + float(intercept))
    
    return err

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
usg = """%prog [options] slope intercept input_file output_file
                                  OR
       %prog [options] -A input_file output_file"""
                
parser=optparse.OptionParser(usage=usg)
parser.add_option("-F","--pjFOV",dest="pjFOV",
                  help="The FOV of the projections")
## No single offset.  In this case, we need to add 
## the intercept and the slope
#parser.add_option("-o","--offset",dest="offset",
#                  help="The offset of the rotational axis from the center of the detectors (in units of pixels")
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
                  help="view antialiasing factor")
parser.add_option("-A","--auto",dest="auto",default="False",action="store_true",
                  help="automatic offset calculation")

(options, args) = parser.parse_args()

verbose = options.verbose

if( (options.auto and len(args) != 2) or (not options.auto and len(args) != 4)):
    parser.error('Not enough arguments.')

if(options.auto):
    infile = args[0]
    outfile = args[1]
else:
    slope = float(args[0])
    intercept = float(args[1])
    infile = args[2]
    outfile = args[3]

if(not os.path.exists(infile)):
    parser.error("The specified input file %s does not exist." % infile)

if(os.path.exists(outfile) and not options.clobber):
    parser.error("The specified output file %s already exists." % outfile)

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

## If we're auto, then time to calculate the slope and the intercept
if(options.auto):
    coarse_offset = get_coarse_offset(infile)
    if(verbose):
        print "Coarse offset is", coarse_offset
    sls = (start_slice +  array([1/5., 2/5., 3/5., 4/5.])  * (end_slice-start_slice)).astype(int16)
    offsets = zeros(shape(sls),float32)
    for i in range(shape(sls)[0]):
        if(verbose):
            print "Now calculating offset for slice", sls[i]
        offsets[i] = float(get_fine_offset(infile,sls[i],0.5,start_offset=coarse_offset-5.0,end_offset=coarse_offset+5.0))
        
    if(verbose):
        print "Output is:", sls, offsets
    ## At this point we should have a sample of the fine offset at each position.
    ## Now we need to fit it to a straight line
    initslope = 1.*(offsets[-1] - offsets[0])/(sls[-1] - sls[0])
    initintercept = 1.*offsets[-1] - 1.*initslope*sls[-1]
    plsq = leastsq(straightline, [ initslope, initintercept ], (offsets, sls))
    slope = plsq[0][0]
    intercept = plsq[0][1]
    if(verbose):
        print "Calculated slope=",slope,"intercept=",intercept
        
if(verbose):
    print "Reconstructing with slope %2.4f and %2.2f" % (slope, intercept)

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

if(verbose):
    print "Minc dimensions are", views, slices, dets
    
td = make_tempdir()
if(verbose):
    print "Working in", td

## Necessary only because of the silliness of python, see the for loop for documentation
#offset_str = "-o %f" % offset
pjFOV_str = "-F %f" % pjFOV
offsetview_str = "-v %d" % offsetview
antialias_str = "--a %s" % options.antialias

if(ROI):
    ROI_str = "-r %s" % ROI
else:
    ROI_str = ""

if(reconsize):
    reconsize_str = "-r %s" % options.reconsize
else:
    ## figure out maximum reconsize
    maxdets = int(ceil(dets+2*max(abs(start_slice*slope+intercept),abs(end_slice*slope+intercept))))
    reconsize_str = "-r %s,%s" % (maxdets, maxdets)

for i in range(start_slice, end_slice):
    execute_command("mincreshape -clobber %s -start 0,%d,0 -count %d,1,%d %s %s/slice%04d.mnc" % 
                           (get_verbose_minc_option(verbose), i, views, dets, infile, td, i))
    ## Can't do this because python does not release memory, and so I run out of memory
    ## See: http://evanjones.ca/python-memory.html -- also, note that this is fixed in python 2.5
    #reconstruct("%s/slice%04d.mnc" % (td, i), "%s/recon%04d.mnc" % (td, i), offset, reconsize, pjFOV, ROI, offsetview)
    offset_str = "-o %2.2f" % (i*slope + intercept) 
    execute_command("python /home/jwalls/Australiawork/packages/opt_tools/individrecon.py %s %s %s %s %s %s %s/slice%04d.mnc %s/recon%04d.mnc" %
                    (antialias_str, offset_str, pjFOV_str, offsetview_str, reconsize_str, ROI_str, td, i, td, i))
    

execute_command("mincconcat -clobber %s -concat_dimension zspace -sequential %s/recon* %s" % 
                       (get_verbose_minc_option(verbose), td, outfile))

cleanup(td)
