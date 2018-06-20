
from scipy import *
import py_minc
from scipy.optimize import leastsq

def gaussian(parameters, y, x):

	a, b, c, d, e, f = parameters
	asdf = e* exp(-((x-a)/b)**2)*exp(-((x-c)/d)**2) - f
	err = y -  asdf
	return err


filen = 'exp7/sl0.mnc'
a = py_minc.ArrayVolume(filen,typecode=int16)
views,null,dets = a.get_sizes()
im1 = squeeze(a.array)[700:800,100:200]

filen = 'exp7/sl399.mnc'
a = py_minc.ArrayVolume(filen,typecode=int16)
views,null,dets = a.get_sizes()
im2 = squeeze(a.array)[700:800,100:200]

dlist = ravel(im1).tolist()
mx = max(ravel(im1))
flr = mean(ravel(im1))
ind = dlist.index(mx)

indx = ind/100
indy = mod(ind,100)

plsq = leastsq(gaussian, [indx,10,indy,10,mx,flr], (im1, arange(-50,50), arange(-50,50)))


## indices = []
## widths = []
## amps = []

## for i in range(views):
## 	dlist = data[i,:].tolist()
## 	mx = max(data[i,:])
## 	flr = mean(data[i,:])
## 	ind = dlist.index(mx)
## 	if(i>200):
## 		ind = 951 - int(floor(indices[i-200]))
## 	print ind, mx, flr
## 	plsq = leastsq(gaussian, [ind, 10, mx, flr], (data[i,ind-50:ind+50], arange(ind-50,ind+50)))
## 	print i, plsq[0]
## 	indices.append(plsq[0][0])
## 	widths.append(plsq[0][1])
## 	amps.append(plsq[0][2])

## print "Done!"
