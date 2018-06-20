
## All this file does is set the variables for use in the fdprecon.py file

## CHANGE THESE VARIABLES FOR YOUR APPLICATION

#data = "/projects/mice/jwalls/data/opt/fdp_simulation/actualdata/test/data.mnc"
#output = "/projects/mice/jwalls/data/opt/fdp_simulation/actualdata/test/recon.mnc"

data = "/projects/mice/jwalls/data/opt/mice_data/fluobeads/MICe01_438_1_GFP1_gbb1/data.mnc"

FOV = 3.3411

## These variables are fixed - no need to change

import tempfile
import os

workdir = tempfile.mktemp('-fdprecon-%s' % os.getlogin())
phantom1 = "%s/phantom1.phm" % workdir
proj1 = "%s/projections1.pj" % workdir
phantom2 = "%s/phantom2.phm" % workdir
proj2 = "%s/projections2.pj" % workdir
variance_iffile = "%s/outfile.if" % workdir

PHM2PJ = "/micehome/workspace/CTSim/tools/ctsimtext phm2pj"
PJREC =  "/micehome/workspace/CTSim/tools/ctsimtext pjrec"

offsets = range(-80,81)
stds = []
subpix_stepsize = 0.1

## THESE ARE ONLY TEMPORARY TEST VARIABLES, DON'T USE THESE IF REALLY RUNNING SCRIPT

views = 400
dets = 1036
slices = 1360
reconsize = 1036

from scipy import *

ff = zeros((views,dets),float32)

from tomography import projections
pj = projections.Projections()
pj.read("/tmp/mytest.pj")
dat = pj.get_data()
g1 = shape(dat)[1]
g2 = dets

left = (g1 - g2)/2  ## Note: These will always be even
right = (g1 + g2)/2 ## as it is either odd-odd or even-even

ff = dat[:,left:right]

import py_minc

a = py_minc.ArrayVolume("exp7/slice746.mnc")

ff = squeeze(a.array).astype(float32)
