
# testing of FDR filtering
import os
from numpy import *
import scipy
from scipy.misc import *
import pylab
from scipy import fftpack
from pyminc.volumes.factory import *

import fdr_correction

# psf coordinate system
"""
#***** try simulating psf

sizes = (1024, 785)
steps = (0.0127, 0.0127)
starts = (-512*steps[0], -384.71*steps[1])
centre = (-0.5, -0.5)
width = 0.01
depth_dependent_width = 0.02

psf = fdr_correction.simulate_2d_gaussian_psf(sizes, steps, starts, centre, width, depth_dependent_width)
pylab.imshow(psf.array[0], extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot)
pylab.xlabel('distance from detector')
pylab.ylabel('detector position')
pylab.title('simulated psf')
pylab.show()

#*********
"""
infile = volumeFromFile('psf512_order.mnc')


psf = infile.data
psf_sizes = (infile.sizes[0],infile.sizes[1],infile.sizes[2])
psf_steps = (infile.data.separations[0],infile.data.separations[1],infile.data.separations[2])
infile.closeVolume()


#psf_starts = (0, -psf_sizes[1]*psf_steps[1]/2.0, -psf_sizes[2]*psf_steps[2]/2.0)
# used with 'new_psf_511_order.mnc' psf_starts = (0, -psf_sizes[1]*psf_steps[1]/2.0, -377.5*psf_steps[2]-5)
#psf_starts = (0, -511*psf_steps[1], -377.5*psf_steps[2]-5)
# read from Display using higest intensity
#psf_starts = (0, -511.21*psf_steps[1],-384.71*psf_steps[2]+(217.80/2.0))
#psf_starts = (0, -511.21*psf_steps[1],-384.71*psf_steps[2])
psf_steps=(0.0127,0.0127,0.0127)

psf_starts = (0, -512*psf_steps[1],-384.71*psf_steps[2])
psf_sizes=(1,1024,785)

# used with crop.mnc psf_starts = (0, -psf_sizes[1]*psf_steps[1]/2.0, -psf_sizes[2]*psf_steps[2]/2.0+30)



npsf = psf.array
#npsf = psf
#***********

"""
psf_sizes = (1, 1024, 785)
psf_steps = (2.45, 2.45, 2.45)
#psf_starts = (-2506.35,-1920.8)
#centre = (0.0, 0.0)
centre = (0.0, 0.0)
psf_starts = (0,-1254.4,-961.625)
width = 5.0
depth_dependent_width = 0.05



twod_psf_sizes =(psf_sizes[1], psf_sizes[2])
twod_psf_starts =(psf_starts[1],psf_starts[2])
twod_psf_steps =(psf_steps[1],psf_steps[2])

psf = fdr_correction.simulate_2d_gaussian_psf(twod_psf_sizes, twod_psf_steps, twod_psf_starts, centre, width, depth_dependent_width)

outfile = volumeFromDescription('simulated_psf.mnc',('xspace', 'yspace', 'zspace'), (1, psf_sizes[1], psf_sizes[2]),(0,psf_starts[1],psf_starts[2]),(1, psf_steps[1],psf_steps[2]),dtype='float',volumeType='short')
outfile.data=psf.array
outfile.writeFile()
outfile.closeVolume()

npsf = psf.array

"""

print psf_starts,psf_steps,psf_sizes

zero_offset = npsf.shape[1]/2

fig1=pylab.figure()
pylab.imshow(npsf[0], extent = (psf_starts[2], psf_starts[2]+(psf_sizes[2]-1)*psf_steps[2],psf_starts[1], psf_starts[1]+(psf_sizes[1]-1)*psf_steps[1]), cmap = pylab.cm.hot, aspect=1)
pylab.xlabel('distance from detector')
pylab.ylabel('detector position')
pylab.title('psf')
pylab.savefig('psf.pdf')
#pylab.show()

new_psf_to_write= npsf.copy()
#new_psf_to_write.shape=(1,psf_sizes[1],psf_sizes[2])
outfile = volumeFromDescription('psf3.mnc',('zspace', 'yspace', 'xspace'), (1,psf_sizes[1], psf_sizes[2]),(0,psf_starts[1],psf_starts[2]),(1, psf_steps[1],psf_steps[2]),volumeType='float',dtype='float') 
outfile.data = new_psf_to_write
outfile.writeFile()
outfile.closeVolume()


# sinogram coordinates (setting sinogram start, step, size as 2D)
# note : sizes[1] # of projection angles --> 2pi/sizes[1]

infile = volumeFromFile('slice512_order.mnc')
sino = infile.data[0]
sizes = (infile.sizes[1],infile.sizes[2])
sizes_final = sizes
steps = (infile.data.separations[1],2*pi/sizes[1])
steps_final = steps
starts = (infile.data.start[1],0.0)
starts_final = starts
infile.closeVolume()

"""
# flipping sinogram for Display 
fliped_sino = flipud(fliplr(sino))
print fliped_sino.shape
fig2=pylab.figure()
pylab.imshow(fliped_sino , extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=0.005*1.5)
pylab.xlabel('projection angle (radians)')
pylab.ylabel('detector position')
pylab.title('sinogram')
#pylab.show()
pylab.savefig('sinogram.pdf')

# saving sinogram for back projection in matlab
savetxt('sinogram.txt', fliped_sino, delimiter=',')
b = fliped_sino.copy()
b.shape=(1,sizes[0], sizes[1])
outfile = volumeFromDescription('sin.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float') 
outfile.data = b
outfile.writeFile()
outfile.closeVolume() 


def fftshift(a):
    return roll(roll(a, a.shape[0]/2, axis = 0), a.shape[1]/2, axis=1) 

# coordinates for fft of sinogram

FOV = (1/float(steps[0]), 1/float(steps[1]))
starts = (-FOV[0]/2.,-FOV[1]/2.)
steps = (FOV[0]/float(sizes[0]),FOV[1]/float(sizes[1]))


# show fourier transform of sinogram
f_sino = fftpack.fft2(sino)
new_f_sino = fftshift(abs(f_sino))


fig3=pylab.figure()
pylab.imshow(new_f_sino, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot,aspect = 200)
pylab.xlabel('angular frequency')
pylab.ylabel('detector spatial frequency')
pylab.title('magnitude of 2D fft of sinogram')
pylab.savefig('fft_sino.pdf')

new_f_sino_to_write= new_f_sino.copy()
new_f_sino_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('fft_2d.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float') 
outfile.data = new_f_sino_to_write
outfile.writeFile()
outfile.closeVolume() 

print npsf.shape
new = roll(fftpack.fft(npsf[0], axis=0) , npsf.shape[1]/2, axis=0)
more_new=abs(new)
#print 'more_new',more_new.shape
more_new.shape=(1,1024,785)
outfile = volumeFromDescription('psf_fft_shift.mnc',('zspace', 'yspace', 'xspace'), (1,1024, 785),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float') 
outfile.data = more_new
outfile.writeFile()
outfile.closeVolume() 

old = fftpack.fft(npsf[0], axis=0)
more_old=abs(old)
more_old.shape=(1,1024,785)
outfile = volumeFromDescription('psf_fft.mnc',('zspace', 'yspace', 'xspace'), (1,1024, 785),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float') 
outfile.data = more_old
outfile.writeFile()
outfile.closeVolume() 

"""

# coordinates for FDR filter

h = fdr_correction.psf_in_fdr_space(sizes, steps, starts, npsf, psf_starts,psf_steps,psf_sizes)

new_h = abs(h)
"""

fig4=pylab.figure()
pylab.imshow(new_h, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect =200)
pylab.xlabel('angular frequency')
pylab.ylabel('detector spatial frequency')
pylab.title('psf transformed to fdr space')
pylab.savefig('fdr.pdf')


new_h_to_write = new_h.copy()
new_h_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('fdr.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = new_h_to_write
outfile.writeFile()
outfile.closeVolume()

#wiener filter (meant to get rid of noise and avoid division by zero
# or very small values which create super large results. changed constant from 0.05-->0.5
# result is much cleaner
n = max(abs(h).flat)*0.5
print max(abs(h).flat)
# n = 3259.01961179
# n**2 = 10621208.830031842
# conj is the complex conjugate
# The complex conjugate of a complex number is obtained 
# by changing the sign of its imaginary part. 1+2j --> 1-2j
# creating inverse of fdr filter
hinv = h.conj()/(abs(h)**2 + n**2)

new_hinv = abs(hinv)


fig5=pylab.figure()
pylab.imshow(new_hinv, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=200)
pylab.xlabel('angular frequency')
pylab.ylabel('detector spatial frequency')
pylab.title('deconvolution filter in fdr space')
pylab.savefig('deconv_fdr.pdf')
#pylab.show()



new_hinv_to_write = new_hinv.copy()
new_hinv_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('fdr_inverse.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = new_hinv_to_write
outfile.writeFile()
outfile.closeVolume()


#calculate the limit recovery filter |Hinv_limit|
# coordinates for FDR filter
fdrlimit=5
limit_recovery_filter = fdr_correction.limit_recovery_filter(sizes, steps, starts, hinv, fdrlimit)

new_limit_recovery = abs(limit_recovery_filter)
"""
fig6=pylab.figure()
pylab.imshow(new_limit_recovery, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=200)
pylab.xlabel('angular frequency')
pylab.ylabel('detector spatial frequency')
pylab.title('max limit recovery filter')
pylab.savefig('max_limit_recovery.pdf')
"""

new_limit_recovery_to_write = new_limit_recovery.copy()
new_limit_recovery_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('limit_recovery_filter.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = new_limit_recovery_to_write
outfile.writeFile()
outfile.closeVolume()

"""
#calculate roll-off fitler
weight=100
maxslope=0.95
roll_off_filter = fdr_correction.roll_off_filter(sizes, steps, starts,weight,maxslope)

fig7=pylab.figure()
pylab.imshow(roll_off_filter, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=200)
pylab.xlabel('')
pylab.ylabel('')
pylab.title('')
pylab.savefig('roll_off_filter.pdf')

roll_off_to_write = roll_off_filter.copy()
roll_off_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('roll_off_filter.mnc',('yspace', 'xspace', 'zspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = roll_off_to_write
outfile.writeFile()
outfile.closeVolume()

#calculate band_limit fitler
bandlimit=0.9
band_limit_filter = fdr_correction.band_limit_filter(sizes, steps, starts,bandlimit)

fig8=pylab.figure()
pylab.imshow(band_limit_filter, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=200)
pylab.xlabel('')
pylab.ylabel('')
pylab.title('')
pylab.savefig('band_limit_filter.pdf')

band_limit_to_write = band_limit_filter.copy()
band_limit_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('band_limit_filter.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = band_limit_to_write
outfile.writeFile()
outfile.closeVolume()
"""
# coordinates for restored sinogram --> sizes_final, steps_final, starts_final

sizes = sizes_final
steps = steps_final
starts = starts_final

#applying fdr filter only
rest_sino = fftshift(fftpack.ifft2(fftshift(f_sino) *h))

new_rest_sino = abs(rest_sino)
"""

fig9=pylab.figure()
pylab.imshow(new_rest_sino, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=0.001*8)
pylab.xlabel('')
pylab.ylabel('')
pylab.title('restored sinogram --> fdr filter only')
pylab.savefig('h_sinogram.pdf')
"""

new_rest_sino_to_write = new_rest_sino.copy()
new_rest_sino_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('h_sinogram.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = new_rest_sino_to_write
outfile.writeFile()
outfile.closeVolume()

#applying inverse fdr filter only
rest_sino = fftshift(fftpack.ifft2(fftshift(f_sino) * hinv ))

new_rest_sino = abs(rest_sino)
"""
fig10=pylab.figure()
pylab.imshow(new_rest_sino, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=0.008)
pylab.xlabel('')
pylab.ylabel('')
pylab.title('restored sinogram --> inverse fdr filter only')
pylab.savefig('hinv_sinogram.pdf')
"""

new_rest_sino_to_write = new_rest_sino.copy()
new_rest_sino_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('hinv_sinogram.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = new_rest_sino_to_write
outfile.writeFile()
outfile.closeVolume()
#pylab.show()

#applying max limit recovery filter

rest_sino = fftshift(fftpack.ifft2(fftshift(f_sino) * limit_recovery_filter ))

new_rest_sino = abs(rest_sino)
"""
fig11=pylab.figure()
pylab.imshow(new_rest_sino, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=0.008)
pylab.xlabel('')
pylab.ylabel('')
pylab.title('restored sinogram --> max recovery filter only')
pylab.savefig('hlimit_sinogram.pdf')
"""

new_rest_sino_to_write = new_rest_sino.copy()
new_rest_sino_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('hlimit_sinogram.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = new_rest_sino_to_write
outfile.writeFile()
outfile.closeVolume()


"""
#applying complete frequency space filter
rest_sino = fftshift(fftpack.ifft2(fftshift(f_sino) * limit_recovery_filter * roll_off_filter * band_limit_filter))

new_rest_sino = abs(rest_sino)

fig12=pylab.figure()
pylab.imshow(new_rest_sino, extent = (starts[1], starts[1]+(sizes[1]-1)*steps[1],starts[0], starts[0]+(sizes[0]-1)*steps[0]), cmap = pylab.cm.hot, aspect=0.008)
pylab.xlabel('')
pylab.ylabel('')
pylab.title('restored sinogram --> complete frequency filter')
pylab.savefig('restored_sinogram.pdf')

pylab.show()

new_rest_sino_to_write = new_rest_sino.copy()
new_rest_sino_to_write.shape=(1,sizes[0],sizes[1])
outfile = volumeFromDescription('restored_sinogram.mnc',('zspace', 'yspace', 'xspace'), (1,sizes[0], sizes[1]),(0,starts[0],starts[1]),(1, steps[0],steps[1]),volumeType='float',dtype='float')
outfile.data = new_rest_sino_to_write
outfile.writeFile()
outfile.closeVolume()

#save as text for back projection algorithm
savetxt('restored_sinogram.txt', new_rest_sino, delimiter=',')
"""
