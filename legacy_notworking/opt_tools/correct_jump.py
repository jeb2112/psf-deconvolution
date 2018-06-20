#########################
# reconall.py
# Takes a multislice minc file projection
# plugs it through reconstruction
#########################

from scipy import *
import tempfile, os, sys, optparse, py_minc
import libOPT

##################################################################################
####### SCRIPT STARTS HERE #######################################################
##################################################################################

usg = 'Usage: %prog [options] input_file output_file'
parser=optparse.OptionParser(usage=usg)
parser.add_option("-j","--jump",dest="jump",default="auto",
          help="The number of pixels (including subpixels) that the sample has jumped.  Use auto for automatic jump determination.")
parser.add_option("-q","--quiet",dest="verbose",default="True",
          action="store_false",help="Don't echo the command")
parser.add_option("-c","--clobber",dest="clobber",
          action="store_true",help="Overwrite the destination file")
parser.add_option("-t","--testonly",dest="testonly",default=False,action="store_true",
                  help="Test only, don't actually output anything.")

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

if(options.jump == "auto"):
  jump = libOPT.get_coarse_jump(infile)
  # jump = get_fine_jump(infile,int((end_slice-start_slice)/2+start_slice),float(options.accuracy),
  #             start_jump=jump_guess-float(options.window)/2.0,
  #             end_jump=jump_guess+float(options.window)/2.0,
  #             plotfile=options.plotfile)
else:
  jump = float(options.jump)

if(options.verbose):
  print "Reconstructing with jump %2.2f" % jump

if(options.testonly):
  sys.exit()

libOPT.correct_jump(infile,outfile,jump)

if(options.verbose):
  print "Done!"
