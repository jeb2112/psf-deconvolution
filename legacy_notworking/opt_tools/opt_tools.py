#!/usr/bin/env python

###################################
#
# opt_tools.py
# python module to handle OPT tasks
#
# Johnathon R. Walls 2006
#
# TODO:
# - lots
###################################

from scipy import *
import py_minc
import tempfile
import os
import shutil

verbose = 1

if(__name__ == "__main__"):
	print "Don't run me"

def execute_command(str):
	if(verbose):
		print "Executing: ", str
	os.system(str)

def tempfile_handler(dir):
	if(os.uname()[0] == 'IRIX64'):
		tempfile.tempdir = '/home/%s' % os.getlogin()
	elif (os.environ.get('TMPDIR')):
		tempfile.tempdir = os.environ.get('TMPDIR')
	else:
		tempfile.tempdir = '/tmp/'
	workdir = tempfile.mktemp('opt_tools')
	if(dir):
		os.mkdir(workdir)
	return workdir


def get_tempfile():
	return tempfile_handler(0)

def make_tempdir():
	return tempfile_handler(1)

def cleanup_tempdir(f):
	shutil.rmtree(f)

def cleanup_tempfile(f):
	os.remove(f)
	
def get_tiff_size(filein):
	
	f1 = get_tempfile()
	cmd = "identify %s > %s" % (filein, f1)
	execute_command(cmd)
	f = open(f1)
	g = f.readlines()[0].split(' ')[2].split('x')
	f.close()
	cleanup_tempfile(f1)
	return g

def tifftominc(filein, fileout):

	im_min = 0
	im_max = 4095
	imsize = get_tiff_size(filein)
	cols = imsize[0]
	rows = imsize[1]

	f = get_tempfile()
	cmd = 'convert -type Grayscale -size %sx%s %s %s' % \
                          (cols, rows, filein, f)
	execute_command(cmd)
	cmd = "cat %s | rawtominc -clobber -short -unsigned -xstep 1 -ystep 1 " \
	      "-zstep 1 -ct -origin 0 0 0 -range %f %f -real_range %f %f %s 1 %s %s" % \
          (f, im_min, im_max, im_min, im_max, fileout, rows, cols)
	execute_command(cmd)
	cleanup_tempfile(f)

