#!/usr/bin/env python

import libOPT
import optparse
import os

## First we move the image files to the views directory

usg = "Usage: %prog [options] directory"
parser = optparse.OptionParser(usage=usg)
parser.add_option("-q","--quiet",dest="verbose",default=False,action="store_false",
                  help="chatty?")

(options, args) = parser.parse_args()
       
if(len(args)!=1):
    parser.error("Incorrect arguments")

d = args[0]

if(not os.path.exists(d)):
    parser.error("The directory %s does not exist" % d)
    
if(not os.path.isdir(d)):
    parser.error("The argument is not a directory")
    
if(not os.path.exists("%s/views" % d)):
    os.mkdir("%s/views" % d)
libOPT.execute_command("mv %s/image_* %s/views" % (d,d))
if(not os.path.exists("%s/binned" % d)):
   os.mkdir("%s/binned" % d)

for i in range(0,400):
    libOPT.tifftominc("%s/views/image_%04d" % (d,i),
                      "%s/views/v%04d.mnc" % (d,i))
    libOPT.execute_command("python %s/rebindata.py -b 2 %s/views/v%04d.mnc %s/views/b_%04d.mnc" % (libOPT.WORKDIR, d,i,d,i))

libOPT.execute_command("mincconcat -sequential -valid_range 0 4095 -float -concat_dimension zspace %s/views/v* %s/data_1.mnc" % (d,d))
libOPT.execute_command("mincconvert -2 %s/data_1.mnc %s/data.mnc" % (d,d))
libOPT.execute_command("mincconcat -sequential -valid_range 0 16380 -float -concat_dimension zspace %s/views/b* %s/binned/data_1.mnc" % (d,d))
libOPT.execute_command("mincconvert -2 %s/binned/data_1.mnc %s/binned/data.mnc" % (d,d))


