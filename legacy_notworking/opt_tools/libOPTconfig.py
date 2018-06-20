#!/usr/bin/env python

import os

## This script doesn't do anything other than just set some variables

if(__name__=="main"):
    print "Nope, can't run this."

## Change below according to your setup
WORKDIR = "/micehome/jwalls/OPT_8.04/opt_tools/"
CTSIMDIR="/micehome/jwalls/OPT_8.04/CTSim_amd64/"
tmpdir = os.getenv("LIBOPT_TMPDIR")
if(tmpdir is None):
    WORKINGDIR = "/tmp/"
else:
    WORKINGDIR = tmpdir
    
## You can change these if you want too.
PHM2PJ = "%s/tools/ctsimtext phm2pj" % CTSIMDIR
PJREC = "%s/tools/ctsimtext pjrec" % CTSIMDIR
