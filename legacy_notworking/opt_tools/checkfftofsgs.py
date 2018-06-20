
from scipy import *
import py_minc


fls = [ "exp2/slice448.mnc", "exp3/slice293.mnc", "exp4/slice286.mnc", "exp5/slice299.mnc",
		"exp6/slice1077.mnc",  "exp7/slice746.mnc",  "exp8/slice752.mnc", "exp9/slice793.mnc" ]

psf_ffts = zeros((8,400,1036),Complex64)
sgs = zeros((8,400,1036),float32)

i = 0

for fl in fls:
	print "Analyzing", fl
	a = py_minc.ArrayVolume(fl,typecode=float32)
	b = squeeze(a.array)
	sgs[i,:,:] = b
	for j in range(400):
		psf_ffts[i,j,:] = fftshift(fft(fftshift(b[j,:])))
	i = i + 1

inv = 1

def showwhere(im,fact=1.0):
	if(inv):
		imshow(-where(im>fact*max(ravel(im)),fact*max(ravel(im)),im))
	else:
		imshow(where(im>fact*max(ravel(im)),fact*max(ravel(im)),im))

def showmag(num, fact=1.0):
	im = abs(squeeze(psf_ffts[num-2,:,:]))
	showwhere(im,fact)

def showabs(num,fact=1.0):
	showmag(num,fact)

def showphase(num, fact=1.0):
	im = angle(squeeze(psf_ffts[num-2,:,:]))
	showwhere(im,fact)

def showreal(num,fact=1.0):
	im = real(squeeze(psf_ffts[num-2,:,:]))
	showwhere(im,fact)

def showimag(num,fact=1.0):
	im = imag(squeeze(psf_ffts[num-2,:,:]))
	showwhere(im,fact)

def showsg(num,fact=1.0):
	im = squeeze(sgs[num-2,:,:])
	showwhere(im,fact)
