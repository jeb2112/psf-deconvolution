from scipy import *
import os

im_min = 0
im_max = 4095
datadir = '.'
outfile = 'data.mnc'
filelen = 400
rows = 1360
cols = 1036


for i in range(400):
	os.system("convert -type Grayscale -size 1036x1360 image%04d image%04d.gray" %  (i,i))


if(os.uname()[0] == "IRIX64"):
	datadir = os.getcwd()
	cmd = "cat *.gray | rawtominc -short -unsigned -xstep 1 -ystep 1 -zstep 1 -ct -origin 0 0 0 -range %f %f -real_range %f %f %s/%s %s %s %s" % \
		  (im_min, im_max, im_min, im_max, datadir, outfile, filelen, rows, cols)
	os.system(cmd)
	
else:
	cmd = "cat *.gray | dd conv=swab | rawtominc -short -unsigned -xstep 1 -ystep 1 -zstep 1 -ct -origin 0 0 0 -range %f %f -real_range %f %f %s %s %s %s" % \
		  (im_min, im_max, im_min, im_max, outfile, filelen, rows, cols)
	os.system(cmd)
	print "Command executed:", cmd


		
