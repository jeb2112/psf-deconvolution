
from scipy import *
import py_minc
from scipy.optimize import leastsq

import sys

def gaussian(parameters, y, x):

	a, b, c, d = parameters
	err = y -  float(c)* exp(-((x-float(a))/float(b))**2) - float(d)

	return err


filen = sys.argv[1]
a = py_minc.ArrayVolume(filen,typecode=float32)
views,dets = shape(squeeze(a.array))
data = squeeze(a.array)

#myhist = Histogram(data.flat, 500)
#ind = compress(equal(myhist[:,1],max(myhist[:,1])),myhist[:,0])
## Now get the fitted value
#fit = leastSquaresFit(gaussian, (ind[0], 10, max(myhist[:,1])), myhist.array)
## fit[0][0] is the gaussian peak
#print fit[0][0]
## At this point we could check the goodness, try again with a new initial
## guess, or just give up if guess is no good

indices = []
widths = []
amps = []

for i in range(views):
	dlist = data[i,:].tolist()
	mx = max(data[i,:])
	flr = mean(data[i,:])
	ind = dlist.index(mx)
	mx = mx - flr
	#if(i>200):
	#	ind = 951 - int(floor(indices[i-200]))
	print ind, mx, flr
	plsq = leastsq(gaussian, [float(ind), 10., float(mx), float(flr)], (data[i,ind-20:ind+20], arange(ind-20,ind+20)))
	print i, plsq[0]
	indices.append(plsq[0][0])
	widths.append(plsq[0][1])
	amps.append(plsq[0][2])

print "Done!"

def loop2(arr):
	sz = shape(arr)
	ret = zeros((sz[0]*2,),float32)
	ret[0:sz[0]] = asarray(arr).astype(float32)
	ret[sz[0]:] = asarray(arr).astype(float32)
	return ret

gplt.plot(widths)
gplt.title(filen)

def check_fit(num):

	i = num

	dlist = data[i,:].tolist()
	mx = max(data[i,:])
	flr = mean(data[i,:])
	ind = dlist.index(mx)

	mx = mx-flr
	gg = mx * exp(-((arange(1036)-ind)/10.)**2) + flr
	plsq = leastsq(gaussian, [float(ind), 10., float(mx), float(flr)], (data[i,ind-50:ind+50], arange(ind-50,ind+50)))
	ind2 = plsq[0][0]
	wid2 = plsq[0][1]
	mx2 = plsq[0][2]
	flr2 = plsq[0][3]
	gg2 = mx2* exp(-((arange(1036)-ind2)/wid2)**2) + flr2
	print ind, 10., mx, flr
	print plsq[0]

	gplt.plot(arange(1036),data[i,:],arange(1036),gg,arange(1036),gg2)
