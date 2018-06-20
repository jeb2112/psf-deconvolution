#!/usr/bin/env python

##############################
#
# prepare_opt_data.py
#
# Johnathon R. Walls, 2006
#
# TODO:
# - brightfield
#
##############################

from libOPT import *
from pylab import *
import os
from optparse import OptionParser
import shutil

#if(__name__="__main__"):
usg = \
'''prepare_opt_data.py can do a few things:
  - prepare brightfield OR fluorescent data
  - process the background value and do the appropriate
    subtraction from the data, depending on method of imaging
  - take into account a single tif stack vs. multiple images
  
To do this it needs to make a few assumptions:
  - All input files are 16 bit files taken with a 12 bit
    camera.  
  - The background file is a single stack
  - averaged files are taken all at the same view rather than
    from subsequent rotations (i.e. im0,im1,im2 are all the
    same view, rather than im0, im400, im800).
  
Good luck!

many files usage: %prog [options] input_dir output_file
stackfile usage:  %prog [options] stackfile outputfile
'''
parser = OptionParser(usage=usg)

parser.add_option("-m", "--method", dest="method",default="fluo",
				  help="fluo or bf for fluorescence and brightfield data.  Default: fluo")
parser.add_option("-b","--bgfile",dest="bgfile",default=False,
          help="The tif file of the background.  -b and -B flags cannot both be active.  Default: False")
parser.add_option("-B","--bgval",dest="bgval",default=False,
          help="The avg background value.  -b and -B flags cannot be both used.  Default: False")
parser.add_option("-a","--average",dest="average",default="1",
          help="How many frames to be averaged. Default: 1")
parser.add_option("-f","--flip",dest="flip",default=False,action="store_true",
          help="rotate file by 90 degrees? Default: False")
parser.add_option("-p","--prefix",dest="prefix",default="image_",help="File prefix.  Default: image_")
parser.add_option("-n","--numviews",dest="numviews",default="400",
          help="Number of views, ignoring averaging.  3 avg or 1 avg, it's still the same number of views.  Default: 400")
parser.add_option("-s","--suffix",dest="suffix",default="",help="File suffix")
parser.add_option("-d","--dirty",dest="clean",action="store_false",default=True,
				  help="Leave behind temp files.  Default: False")
parser.add_option("-q","--quiet",dest="verbose",action="store_false",default=True,
          help="Do not print messages to stdout.  Default: False")
parser.add_option("-c","--clobber",dest="clobber",action="store_true",default=False,
          help="Overwrite existing files.  Default: False")
parser.add_option("-r","--range",dest="rng",default="0,4095",
          help="Range of data..  Default: 0,4095")
parser.add_option("-l","--littleendian",dest="hack_endian",default=False,action="store_true",
          help="Hack for tiffsplit bug.  Default: False")

(options, args) = parser.parse_args()

if(len(args) != 2):
  parser.error("Incorrect arguments.  See help file")

verbose = options.verbose

inputarg = args[0]
outfile = args[1]
outfile = os.path.realpath(outfile)

if(int(options.numviews) <=0):
	parser.error("Error: the number of views specified is invalid")

if(os.path.exists(outfile) and not options.clobber):
  parser.error("Error: the specified output file already exists")

if(options.bgfile and not os.path.exists(options.bgfile)):
  parser.error("Error: the bg file does not exist")

if(options.method != "fluo" and options.method != "bf"):
  parser.error("Incorrect method.  Must be fluo or bf.")
  
if(options.method == "bg" and (not options.bgfile or not options.bgval) ):
  parser.error("You have indicated that this is brightfield data but have not supplied the background file or value")

if(options.bgfile and options.bgval):
  parser.error("Both bgfile and bgval have been specified.  I don't know which one to choose.  Can only use one.")
  
if(options.flip):
  flipstr = "-f"
else:
  flipstr = ""

if(len(options.rng.split(",")) != 2):
  parser.error("The range is incorrectly formatted.")

lower = int(options.rng.split(",")[0])
upper = int(options.rng.split(",")[1])

if(lower > upper):
  parser.error("Range is incorrectly formatted.")

td = make_tempdir()

if(os.path.isfile(inputarg)):
  ## We're doing a tiff stack so let's split it
  split_tiff(inputarg, td, options.prefix, options.hack_endian)
  options.suffix=".tif"
  inputarg = td
  ## This gives us a bunch of tiff files in td
## Let's just make sure everything's there
for i in range(0,int(options.numviews)*int(options.average)):
    if(not os.path.exists("%s/%s%04d%s" % (inputarg, options.prefix, i, options.suffix))):
      parser.error("Error: File %s%04d%s does not exist in the input directory %s" %
            (options.prefix, i, options.suffix, td))  

## Now we can proceed to turn these into minc files,
## whether they are in td or in some other dir
for i in range(0,int(options.numviews)*int(options.average)):
    execute_command("python %s/tifftominc.py %s -r %s %s/%s%04d%s %s/%s%04d.mnc" % 
       (WORKDIR, flipstr, options.rng, inputarg, options.prefix, i, options.suffix, td, options.prefix, i)) 
  
  ## At the end of this loop we have n * A minc files stored in the tempdir td

## Now we average the minc files if necessary

if(int(options.average) > 1):
  for i in range(int(options.numviews)):
    filestoavg = ""
    for j in range(int(options.average)):
      filestoavg = filestoavg + " %s/%s%04d.mnc" % (td, options.prefix, i*int(options.average)+j)
    execute_command("mincmath -add %s %s/tmp%04d.mnc" % (filestoavg, td, i))
    execute_command("mincmath -float -div %s/tmp%04d.mnc -const %s %s/avg%04d.mnc" % (td, i, options.average, td, i))
  options.prefix = "avg"
  
## We now have averaged (if necessary) output files

if(options.bgval):
    # We need to do background subtraction
  if(options.method=="fluo"):
    # We're doing fluorescence
    for i in range(int(options.numviews)):
      execute_command("mincmath -sub %s/%s%04d.mnc -const %5.10f %s/prep%04d.mnc" % 
                     (td, options.prefix, i, float(options.bgval), td, i))
  else:
    # We're doing brightfield, which means we need to do the data transform as well
    for i in range(int(options.numviews)):
      execute_command("minccalc -expression \"log(%5.10f/A[0])\" %s/%s%04d.mnc %s/prep%04d.mnc" % 
      (float(options.bgval), td, options.prefix, i, td, i))
  options.prefix = "prep"

elif(options.bgfile):
  ## First convert the bgfile to a single minc file
  tiffstack_to_minc(options.bgfile, "%s/bg.mnc" % td, flip=options.flip, rang=(lower,upper),hack_endian=options.hack_endian)
  execute_command("mincaverage -avgdim zspace -float -2 %s/bg.mnc %s/avgbg.mnc" % (td, td))
  execute_command("python %s/addzdim.py %s/avgbg.mnc" % (WORKDIR, td))
  
  # We need to do background subtraction
  if(options.method=="fluo"):
    # We're doing fluorescence
    for i in range(int(options.numviews)):
      execute_command("mincmath -sub %s/%s%04d.mnc %s/avgbg.mnc %s/prep%04d.mnc" % (td, options.prefix, i, td, i))
  else:
    # We're doing brightfield, which means we need to do the data transform as well
    for i in range(int(options.numviews)):
      execute_command("minccalc -expression log(A[0]/A[1]) %s/%s%04d.mnc %s/avgbg.mnc %s/prep%04d.mnc" % 
      (td, options.prefix, i, td, td, i))
  options.prefix = "prep"
  
## That's done, now let's concatenate it

execute_command("mincconcat %s -sequential -2 -clobber -concat_dimension zspace %s/%s????.mnc %s" %
          (get_verbose_minc_option(options.verbose), td, options.prefix, outfile))

if(options.clean):
	cleanup(td)
	
