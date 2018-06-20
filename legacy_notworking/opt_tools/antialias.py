from scipy import *
import py_minc
from scipy.weave import converters


def weave_interpmath(oldr,oldc,newr,newc,rows,cols,mat,newmat):
                                                                                                    
	code = """
	float indr,indc,fracr,fracc;
	int findr,findc,i,j;
	for(i=0;i<newr;i++) {
	    for(j=0;j<newc;j++) {
		    indr = rows(i);
			indc = cols(j);
			findr = int(floor(indr));
			findc = int(floor(indc));
			fracr = indr - findr;
			fracc = indc - findc;
			newmat(i,j) = (1-fracr) * (1-fracc) * mat(findr,findc);
			if(findc+1 <= oldc-1)
			    newmat(i,j) += (1-fracr) * (fracc) * mat(findr,findc+1);
                                                                                                    
            if(findr+1 <= oldr-1)
			    newmat(i,j) += (fracr) * (1-fracc) * mat(findr+1,findc);

            if((findr+1 <= oldr-1) && (findc+1 <= oldc-1))
			    newmat(i,j) +=  (fracr)   * (fracc)   * mat(findr+1,findc+1);
        }
    }
                                                                                                    
"""
	weave.inline(code,['oldr','oldc','newr','newc','rows','cols','mat','newmat'],\
				 compiler='gcc',type_converters=converters.blitz)
                                                                                                    


def interp2d(mat,sz):
                                                                                                    
	## Note size can be a multiplier, or a tuple of row/col
	
	if(len(shape(mat)) != 2):
		print "I can only do 2d matrices"
		return mat
	
	if(size(sz) == 2):
		newr = sz[0]
		newc = sz[1]
	else:
		newr = floor(shape(mat)[0] * sz)
		newc = floor(shape(mat)[1] * sz)

        #print newr,newc
	newmat = zeros((newr,newc),mat.typecode())
	ratr = (shape(mat)[0]-1) / (newr-1.)
	ratc = (shape(mat)[1]-1) / (newc-1.)
	#print ratr,ratc
	rows = arrayrange(newr)*(ratr*1.)
	cols = arrayrange(newc)*(ratc*1.)
	#print shape(mat), max(rows), max(cols)
	weave_interpmath(shape(mat)[0],shape(mat)[1],newr,newc,rows,cols,mat,newmat)
	
	return newmat


def build_antialias_mask(mat,viewmax,detmax,view2max,det2max):

	views,dets = shape(mat)
	mask = ones(shape(mat),float32)
	viewstart = 1/pi
	dview = -2./pi/newviews
	detstart = -1/2.
	ddet = 1./newdets

	#print "Starting with", viewmax, detmax, view2max, det2max
	for view in range(views):
		viewspace = viewstart + (view+0.5) * dview
	    #print viewstart, (viewnum + 0.5), dview, (viewnum + 0.5) * dview, viewspace

	for det in range(dets):
		detspace = detstart +(det+0.5) * ddet

	for view in range(views):
		viewspace = viewstart + (view+0.5) * dview
		if(abs(viewspace) > viewmax):
			# roll off with a cosine!
			mask[view,:] = 0
			continue

		for det in range(dets):
			if(abs(detspace) > detmax):
				mask[:,det] = 0
				continue

			detspace = detstart +(det+0.5) * ddet
			#print "coords are", viewspace, detspace
			if(viewspace>view2max and detspace < det2max):
				#print "in R1"
				mask[view,det] = 0
			if(viewspace<-view2max and detspace> -det2max):
				mask[view,det] = 0
				#print "in R2"

	return mask
				

###################################
### MAIN STARTS HERE ##############
###################################

slope = 0.6
max_intercept = 0.3*slope
viewmax = 1/pi/2.
detmax = 0.4
file = "exp3/slice293.mnc"

a = py_minc.ArrayVolume(file)
sg = squeeze(a.array)
sg_fft = fftshift(fft2(fftshift(sg)))

views, dets = shape(sg_fft)

newsg_fft = zeros((views*3,dets),Complex64)
newsg_fft[0:views,:] = sg_fft
newsg_fft[views:views*2,:] = sg_fft
newsg_fft[views*2:views*3,:] = sg_fft

viewsgfft = where(abs(newsg_fft)>7000,7000,abs(newsg_fft))

newviews,newdets = shape(newsg_fft)

mask = zeros(shape(newsg_fft),Int8)
slopes = zeros(shape(newsg_fft),float32)
intercepts = zeros(shape(newsg_fft),float32)
viewstart = 1/pi
dview = -2./pi/newviews
detstart = -1/2.
ddet = 1./newdets

for view in range(newviews):
	viewspace = viewstart + (view+0.5) * dview
	#print viewstart, (viewnum + 0.5), dview, (viewnum + 0.5) * dview, viewspace
	if(abs(viewspace) > viewmax):
		mask[view,:] = 0
		continue
	for det in range(dets):
		detspace = detstart + (det+0.5) * ddet
		if(abs(detspace) > detmax):
			mask[view,det] = 0
			continue
		detspace = detstart + (det+0.5)*ddet
		#print detspace
		intercept = viewspace - slope * detspace
		intercepts[view,det] = intercept
		if(abs(intercept) < max_intercept):
			mask[view,det] = 1
		else:
			mask[view,det] = 0
## 		if(abs(detspace) < 1e-8):
## 			slopes[view,det] = inf
## 		else:
## 			slopes[view,det] = viewspace/detspace

#imshow(mask * viewsgfft)
back = abs(ifft2(mask * newsg_fft))
#imshow(back)
resizesg = interp2d(sg,(1200,1036))

zz = build_antialias_mask(newsg_fft, 1/pi/2., 0.6, 1/pi/2. - 0.05, 0.2)
