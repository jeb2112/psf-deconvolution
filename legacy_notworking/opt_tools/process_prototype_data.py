#
# process_opt_data.py
# Master script for processing data from prototype or bioptonics
# Johnathon R. Walls
# June 2008
#
# The syntax expected is:
# process_opt_data.py directory
# where cirectory contains the OPT Scan data
# The script automatically determines if it's bioptonics data
# or prototype data, if it's brightfield or fluorescence, creates
# a minc file accordingly.  Data can also be FDR-filtered from 
# this script as well as run through the recon.  This results
# in a lot of options to the program.

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

usg = "Usage: %prog [options] directory"
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
parser.add_option("-s", "--scan_output",dest="scan_outputname",default="rawdata.mnc",
                  help="The name of the raw data .mnc file to be written.  Default:  rawdata.mnc")
parser.add_option("-b", "--bg_output",dest="bg_outputname",default="bg.mnc",
                  help="The name of the background data .mnc file to be written.  Default: bg.mnc.  Use none for no background.")
parser.add_option("-p", "--processed_output",dest="processed_outputname",default="data.mnc",
                  help="The name of the processed data .mnc file to be written.  Default: data.mnc")
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
parser.add_option("-x","--bgpath",dest="bgpath",default="/tcphome/jwalls/data/opt/backgrounds",
                  help="The location in which the background images are stored..")
parser.add_option("-y","--psfpath",dest="psfpath",default="/tcphome/jwalls/data/opt/psfs/",
                  help="The location in which the psfs and FDR filters are stored.")

  
(options, args) = parser.parse_args()
if(len(args) != 1):
    parser.error("Incorrect arguments")

outfiles = (options.scan_outputname, options.bg_outputname, options.processed_outputname, options.recon_outputname, 
            options.offsetplot_outputname, options.info_outputname)
fdroutfiles = (options.fdrscan_outputname, options.fdrrecon_outputname, options.fdroffsetplot_outputname)
       
libOPT.verbose = options.verbose

d = os.path.realpath(args[0])

if(not os.path.exists(d)):
  die("The directory %s does not exist." % d)
  
if(not os.path.isdir(d)):
  die("The directory %s is not a directory, skipping ..." % d)

# We're going to do build a script that does the following, in order
# 1. Load up the annotations
# 2. Make sure all the files are there, including background files
# 3. Check to make sure that we have clobber or that the outputs don't already exist
# 4. Convert the scan file to mnc
# 5. If no bg .mnc file, find the tiff, turn it into the appropriate background file
# 6. Link the appropriate background file
# 7. Process the background and data files for final result
# 8. If FDR-filtering is requested, do so
# 9. recon the bg processed data
# 10. recon the fdr-filtered data
# 11. clean everything up (if all of these were successful).

commands_we_ran = []
files_to_remove = []

# 1. Load up the annotations
if(not os.path.exists("%s/ANNOTATION_A.txt" % d)):
  die("The ANNOTATION_A.txt file does not exist in directory %s" % d)
else:
  (status, s) = do_or_die("%s/mac2unix.sh %s/ANNOTATION_A.txt" % (libOPT.WORKDIR, d))
  if(status):
    die(s)
  else:
    say("Found ANNOTATION_A.txt, loading...")

ScanInfo = libOPT.parse_annotations("%s/ANNOTATION_A.txt" % d)
if(ScanInfo is False):
  die("Could not parse ANNOTATION_A.txt")

# 2. Make sure all the files are there, including background files
if(os.path.exists("%s/OPT_Scan.tif" % d)):
  say("Found OPT_Scan.tif, proceeding as a low res scan")
  resolution = "lowres"
  scanfile = "%s/OPT_Scan.tif" % d
else:
  resolution = "highres"
  ## Maybe it's a high res scan?
  for i in range(int(ScanInfo["NumberOfAngles"])):
    scanfile = ()
    infile = "%s/image_%04d" % (d,i)
    if(not os.path.exists(infile)):
      die("Cannot find either OPT_Scan.tif or image_%04d.tif in %s, aborting" % (i,d))
  say("Found files image_0000.mnc through to image_%04d.mnc, proceeding as high res" % (int(ScanInfo["NumberOfAngles"])-1))
  scanfile = d

if(options.bg_outputname != "none"):
  if(os.path.exists("%s/OPT_Background.tif" % d) and not os.path.islink("%s/OPT_Background.tif" % d)):
    say("Found OPT_Background.tif, proceeding as tOPT.")
    mode = "transmission"
  else:
    mode = "emission"
    ## Let's try to find the background file
    bgfile = "%s/%s/%s" % (options.bgpath, libOPT.parse_date(ScanInfo["ScanDate"]), libOPT.get_bg_filename(ScanInfo))
    if(os.path.exists(bgfile)):
      if(os.path.exists("%s/OPT_Background.tif" % d)):
        do_or_die("rm %s/OPT_Background.tif" % d)
        say("Found an old link to the bg file, removing it ...")
      do_or_die("ln -s %s %s" % (bgfile, os.path.realpath("%s/OPT_Background.tif" % d)))
      say("Found the bg file and linked it, proceeding as eOPT")
    else:
      die("Could not find OPT_Background.tif or %s, aborting." % bgfile)
    
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
  
do_or_die("python %s/prepare_opt_data.py -d -c %s %s/%s" % 
          (libOPT.workdir, scanfile, d, options.scan_outputname))

# 5. Process the bg file and the data file
if(options.bg_outputname == "none"):
  if(os.path.exists("%s/%s" % (d, options.processed_outputname))):
    os.system("rm %s/%s" % (d,options.processed_outputname))
  do_or_die("ln -s %s/%s %s/%s" % (d,options.scan_outputname,d,options.processed_outputname))  
else:
  do_or_die("python %s/prepare_background.py -c %s/OPT_Background.tif %s/%s" % 
            (libOPT.WORKDIR, d, d, options.bg_outputname))
  # 6. Process the background and data files for final result
  do_or_die("python %s/process_scan.py -c -m %s %s/%s %s/%s %s/%s" % 
            (libOPT.WORKDIR, mode, d, options.scan_outputname, d, options.bg_outputname,d,options.processed_outputname))
  
# 7. recon the bg processed data
tf = libOPT.get_tempfile()
do_or_die("python %s/reconall.py -c -o twostep -p %s/%s  %s/%s %s" % 
          (libOPT.WORKDIR, d, options.offsetplot_outputname, d, options.scan_outputname, tf))
do_or_die("python %s/convert_for_archive.py -c %s %s/%s" % 
          (libOPT.WORKDIR, tf,d,options.recon_outputname))
libOPT.cleanup(tf)
# 8. If FDR-filtering is requested, do so, along with recon
if(options.fdr_filter):
  fdr_options = ""
  do_or_die("python %s/fdr_filter.py -c %s %s/%s %s/%s" % 
            (libOPT.WORKDIR, fdr_options, d, options.scan_outputname, options.fdrscan_outputname))
  tf = libOPT.get_tempfile()
  do_or_die("python %s/reconall.py -c -o twostep -p %s/%s  %s/%s %s" % 
          (libOPT.WORKDIR, d, options.fdroffsetplot_outputname, d, options.fdrscan_outputname, tf))
  do_or_die("python %s/convert_for_archive.py -c %s %s/%s" % 
          (libOPT.WORKDIR, tf,d,options.fdrrecon_outputname))
  libOPT.cleanup(tf)
  
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
