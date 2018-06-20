#!/usr/bin/env python
#
# process_bioptonics_data.py
# Master script for processing data from bioptonics device
# Johnathon R. Walls
# July 2008
#
# The syntax expected is:
# process_bioptonics_data.py logfile
# where logfile is the .log file of the scan you wish to reconstruct.
# The script automatically determines if it's brightfield or
# fluorescence and creates a minc file accordingly.  
# Data can also be FDR-filtered from this script as well as 
# run through the recon. 

import libOPT
import optparse
import os,sys

def say(s):
  libOPT.provide_feedback(s)

def die(s):
  say(s)
  sys.exit()

def record_commands_to_file():
  f = open("%s/%s" % (d, options.info_outputname))
  f.writelines(commands_we_ran)
  f.close()

def do(s):
  commands_we_ran.append(s)
  return libOPT.execute_command(s)
  
def do_or_die(s):
  (status, s2) = do(s)
  if(status):
    record_commands_to_file()
    die(s2)
  return (status, s2)

def cmd_prep(s):
  return s.replace(" ","\ ")

if(__name__ != "__main__"):
  die("This script can only be run via main")

## First we move the image files to the views directory

usg = "Usage: %prog [options] <logfile>"
parser = optparse.OptionParser(usage=usg)
parser.add_option("-q","--quiet",dest="verbose",default=True,action="store_false",
                  help="Suppress output.  Default: False")
parser.add_option("-c","--clobber",dest="clobber",default=False,action="store_true",
                  help="Overwrite existing files.  Default: False")
parser.add_option("-d","--dirty",dest="dirty",default=False,action="store_true",
                  help="Leave behind original and temporary files")
parser.add_option("-k","--keep_float",dest="archive",default=True,
                  action="store_false",help="Enable this to prevent storing the final recon as short.")
## Options for output file names
parser.add_option("-f","--fdr_filter",dest="fdr_filter",default=False,
                  help="perform FDR-based filtering of the data?  Default: False")
parser.add_option("-s", "--scan_output",dest="scan_outputname",default="data.mnc",
                  help="The name of the raw data .mnc file to be written.  Default:  data.mnc")
parser.add_option("-r", "--recon_output",dest="recon_outputname",default="recon.mnc",
                  help="The name of the recon .mnc file to be written.  Default: recon.mnc")
parser.add_option("-o", "--offsetplot_output",dest="offsetplot_outputname",default="offsetplot.png",
                  help="The name of the auomatically calculated offset plot .png file to be written.  Default: offsetplot.png")
parser.add_option("-i", "--info_output",dest="info_outputname",default="commands.log",
                  help="The name of the commands-executed info text file to be written.  Default: commands.log")
## Options for FDR-filtering output file names
parser.add_option("-S", "--fdrscan_output",dest="fdrscan_outputname",default="fdr_data.mnc",
                  help="The name of the FDR-processed data .mnc file to be written.  Default: fdr_data.mnc")
parser.add_option("-R", "--fdrrecon_output",dest="fdrrecon_outputname",default="fdr_recon.mnc",
                  help="The name of the FDR-processed recon .mnc file to be written.  Default: fdr_recon.mnc")
parser.add_option("-O", "--fdroffsetplot_output",dest="fdroffsetplot_outputname",default="fdroffsetplot.png",
                  help="The name of the FDR-processed automatically calculated offset plot .png file to be written.  Default: fdroffsetplot.png")
## Paths to search
parser.add_option("-y","--psfpath",dest="psfpath",default="/tcphome/jwalls/data/opt/psfs/",
                  help="The location in which the psfs and FDR filters are stored.")

  
(options, args) = parser.parse_args()
if(len(args) != 1):
    parser.error("Incorrect arguments")

outfiles = (options.scan_outputname, options.recon_outputname, 
            options.offsetplot_outputname, options.info_outputname)
fdroutfiles = (options.fdrscan_outputname, options.fdrrecon_outputname, options.fdroffsetplot_outputname)
       
libOPT.verbose = options.verbose

log = os.path.realpath(args[0])

if(not os.path.exists(log)):
  die("The log file %s does not exist." % log)

d = "/".join(log.split("/")[:-1])

if(not os.path.exists(d)):
  die("Could not locate the parent directory %s." % d)
  
# We're going to do build a script that does the following, in order
# 1. Load up the log file
# 2. Make sure all the files are there
# 3. Check to make sure that we have clobber or that the outputs don't already exist
# 4. Convert the scan file to mnc
# 5. recon the raw data
# 6. If FDR-filtering is requested, do so
# 7. recon the fdr-filtered data
# 8. clean everything up (if all of these were successful).

commands_we_ran = []
files_to_remove = []

# 1. Load up the log file
ScanInfo = libOPT.parse_bioptonics_log(log)
if(ScanInfo is False):
  die("Could not parse log file %s" % log)

# 2. Make sure all the files are there

if(ScanInfo["Image Format"] != "TIFF"):
  die("The input file format is not TIFF as expected.")
  
if(ScanInfo["Depth (bits)"] != "16"):
  die("The input file bit depth is not 16 as expected.")

prefix = ScanInfo["Filename Prefix"]
numviews = ScanInfo["Number Of Files"]

infiles = [] 

for i in range(int(numviews)):
  f = "%s/%s%04d.tif" % (d, prefix, i)
  if(not os.path.exists(f)):
    die("Could not find file %s as expected." % f)
  infiles.append(f)

if(ScanInfo["Filter"] == "White(13) light"):
  mode = "transmission"
else:
  mode = "emission"
    
# 3. Check to make sure that we have clobber or that the outputs don't already exist
if(not options.clobber):
  for of in outfiles:
    if(os.path.exists("%s/%s" % (d, of))):
      die("No clobber is active and %s/%s exists, aborting ..." % (d,of))
  if(options.fdr_filter):
    for of in fdroutfiles:
      if(os.path.exists("%s/%s" % (d, of))):
        die("No clobber is active and %s/%s exists, aborting ..." % (d,of))

## We've found the annotation file, the scan file(s) and the background file.  We have room to write
## the files we want.  Time to proceed.  
      
# 4. Convert the scan file to mnc
scanfile = " ".join(infiles)
do_or_die("python %s/tifftominc.py -i -r 0,65535 %s %s/%s" % 
          (libOPT.WORKDIR, scanfile, d, options.scan_outputname))

# 7. recon the raw data
tf = libOPT.get_tempfile()
do_or_die("python %s/reconall.py -c -o twostep -p %s/%s  %s/%s %s" % 
          (libOPT.WORKDIR, d, options.offsetplot_outputname, d, options.scan_outputname, tf))
do_or_die("python %s/convert_for_archive.py -c %s %s/%s" % 
          (libOPT.WORKDIR, tf,d,options.recon_outputname))
libOPT.cleanup(tf)
# 8. If FDR-filtering is requested, do so, along with recon
if(options.fdr_filter):
  say("There is no possible fdr filter for the Bioptonics system.  Sorry!")
#  fdr_options = ""
#  do_or_die("python /micehome/jwalls/workspace/libFDR/fdr_filter.py -c %s %s/%s %s/%s" % 
#            (fdr_options, d, options.scan_outputname, options.fdrscan_outputname))
#  tf = libOPT.get_tempfile()
#  do_or_die("python /micehome/jwalls/workspace/opt_tools/reconall.py -c -o twostep -p %s/%s  %s/%s %s" % 
#          (d, options.fdroffsetplot_outputname, d, options.fdrscan_outputname, tf))
#  do_or_die("python /micehome/jwalls/workspace/opt_tools/convert_for_archive.py -c %s %s/%s" % 
#          (tf,d,options.fdrrecon_outputname))
#  libOPT.cleanup(tf)
  
# 9. clean everything up (if all of these were successful).
if(not options.dirty):
  for file in files_to_remove:
    do_or_die("rm %s" % cmd_prep(file))


#
#
#
#  
#if(not os.path.exists("%s/views" % d)):
#    os.mkdir("%s/views" % d)
#libOPT.execute_command("mv %s/image_* %s/views" % (d,d))
#if(not os.path.exists("%s/binned" % d)):
#   os.mkdir("%s/binned" % d)
#
#for i in range(0,400):
#    libOPT.tifftominc("%s/views/image_%04d" % (d,i),
#                      "%s/views/v%04d.mnc" % (d,i))
#    libOPT.execute_command("python /micehome/jwalls/workspace/opt_tools/rebindata.py -b 2 %s/views/v%04d.mnc %s/views/b_%04d.mnc" % (d,i,d,i))
#
#libOPT.execute_command("mincconcat -sequential -valid_range 0 4095 -float -concat_dimension zspace %s/views/v* %s/data_1.mnc" % (d,d))
#libOPT.execute_command("mincconvert -2 %s/data_1.mnc %s/data.mnc" % (d,d))
#libOPT.execute_command("mincconcat -sequential -valid_range 0 16380 -float -concat_dimension zspace %s/views/b* %s/binned/data_1.mnc" % (d,d))
#libOPT.execute_command("mincconvert -2 %s/binned/data_1.mnc %s/binned/data.mnc" % (d,d))
#
