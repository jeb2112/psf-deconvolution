import libOPT
from optparse import OptionParser
import os
from scipy import *

usg = "Usage: %prog [options] infile"
parser = OptionParser(usage=usg)

parser.add_option("-o","--outfile",dest="outfile",help="Output to a new file",default=False)
parser.add_option("-c","--clobber",dest="clobber",help="Overwrite existing file (applies to --outfile only)", )

(options,args) = parser.parse_args()

if(len(args) != 1):
    parser.error("Incorrect arguments")
    
infile = args[0]

if(not os.path.exists(infile)):
    parser.error("The input file %s does not exist" % infile)

if(options.outfile):
    if(os.path.exists(options.outfile) and not options.clobber):
        parser.error("The specified output file %s already exists." % options.outfile)
    out = os.path.realpath(options.outfile)
else:
    out = os.path.realpath(infile)

a = libOPT.openmincfile(infile)

sz = shape(a.array)
if(len(sz) != 2):
    parser.error("This minc file is improperly formatted -- does not have only 2 dims")
    
b = libOPT.newmincfile((1,sz[0],sz[1]))
b.array[0,:,:] = squeeze(a.array)
b.set_range(min(ravel(a.array)),max(ravel(a.array)))

tf = libOPT.get_tempfile()

b.output(tf)

libOPT.execute_command("mv %s %s" % (tf, out))

