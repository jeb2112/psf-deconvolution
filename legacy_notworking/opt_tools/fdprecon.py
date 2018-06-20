
###########################
# This is the outline
###########################

# 1) read in minc file
# 2) create phantom and projections
# 3) do many reconstruction with various offsets
# 4) measure variance and fit to find maximum
# 5) calculate the appropriate sized phantom and projections
#    - create phantom and projection file
# 6) fit the real data into the new projections
# 7) reconstruct
# 8) read in reconstructions, output to minc file

from variables import *
from scipy import *
import py_minc

def openmincfile(filename, dn=py_minc.FILE_ORDER_DIMENSION_NAMES, tc=int16):
	a = py_minc.ArrayVolume(filename,dim_names=dn, typecode=tc)
	return a

def newmincfile(dims,ord=(py_minc.MIzspace,py_minc.MIyspace,py_minc.MIxspace), \
				ncdt=py_minc.NC_FLOAT,tc=float32):
	a = py_minc.ArrayVolume(dims,ord,nc_data_type=ncdt,typecode=tc)
	return a

def create_phantom(str,file):
	f = open(file,'w')
	f.writelines(str)
	f.close()


def get_pj_data(file):
	from tomography import projections
	pj = projections.Projections()
	pj.read(file)
	a = pj.get_data()
	return a

def write_pj_data(data,file):
	from tomography import projections
	pj = projections.Projections()
	pj.read(file)
	pj.set_data(data)
	pj.write(file)

def get_if_data(file):
	from tomography import imagefile
	iffile = imagefile.ImageFile()
	iffile.fileRead(file)
	a = iffile.get_data()
	return a

def cleanup(d):
	import shutil
	shutil.rmtree(d)

def IsEven(d):
	return (not (d%2))

def execute_command(str):
	print str
	os.system(str)

##################################################################################
####### SCRIPT STARTS HERE #######################################################
##################################################################################

if (not os.path.exists(workdir)):
	os.mkdir(workdir)

## 1) Read in mincfile

#a = openmincfile(data)
#views,slices,dets = a.get_sizes()

detInc = FOV*2/dets
print "detInc=", detInc

dets1 = dets*sqrt(2)
FOV1 = FOV*sqrt(2)

detInc1 = FOV1*2/dets1
print "detInc=",detInc1

dets2 = int(ceil(dets1))
if((IsEven(dets) and not IsEven(dets2)) or (not IsEven(dets) and IsEven(dets2))):
	dets2 = dets2 + 1

FOV2 = FOV1 * dets2/dets1
detInc2 = FOV2*2/dets2
print "detInc=", detInc2
newFOV = FOV * dets2/dets1


## 2) Create phantom and projections

ph1 = "ellipse 0 0 %f %f 0 1" % (FOV, FOV)
create_phantom(ph1,phantom1)
execute_command("%s %s %d %d --rotangle 1 --phmfile %s" % (PHM2PJ, proj1, dets2, views, phantom1))

pjdat = get_pj_data(proj1)

left = (dets2 - dets)/2  ## Note: These will always be even
right = (dets2 + dets)/2 ## as it is either odd-odd or even-even

pjdat[:,left:right] = ff

write_pj_data(pjdat,proj1)

## 3) do many reconstructions with various offsets (pixel offset only)

for offset in offsets:
	execute_command("%s %s %d %d --offset %d" % (PJREC, proj1, variance_iffile, offset, dets, dets, offset))
	ifdata = get_if_data(variance_iffile)
	stds.append(std(ravel(ifdata)))
	print "std is", stds[-1]

## 3a) find most likely range for correct offset

ind = stds.index(max(ravel(stds)))
gplt.plot(stds)

## 4) do many reconstructions with various offsets (subpixel offset only)


offsets2 = arange((offsets[ind]-2),(offsets[ind]+(2+subpix_stepsize)),subpix_stepsize)
stds2 = []

for offset in offsets2:
	execute_command("%s %s %s %d %d --offset %f" % (PJREC, proj1, variance_iffile, dets/2, dets/2, offset))
#	ifdata = get_if_data(variance_iffile)
#	stds2.append(std(ravel(ifdata)))	

## 4a) measure variance and fit to find maximum

tck = interpolate.splrep(offsets2,stds2,s=0)
xnew = arange(offsets2[0],offsets2[-1],0.01)
ynew = (interpolate.splev(xnew,tck,der=0)).tolist()
corroffset = xnew[ynew.index(max(ravel(ynew)))]

# corroffset = offsets2[stds2.index(max(ravel(stds2)))]

print "Correct offset is", corroffset

## 5) calculate the appropriate sized phantom and projections

dts = dets + abs(corroffset*2) + 10

dets1 = (dts)*sqrt(2)
FOV1 = FOV*dts/dets

dets2 = int(ceil(dets1))
if((IsEven(dets) and not IsEven(dets2)) or (not IsEven(dets) and IsEven(dets2))):
	dets2 = dets2 + 1

newFOV = FOV1 * dets2/dets1

## 5a) create phantom and projection file

ph2 = "ellipse 0 0 %f %f 0 1" % (newFOV, newFOV)
create_phantom(ph2,phantom2)
execute_command("%s %s %d %d --rotangle 1 --phmfile %s" % (PHM2PJ, proj2, dets2, views, phantom2))

pjdat = get_pj_data(proj2)

left = (dets2 - dets)/2  ## Note: These will always be even
right = (dets2 + dets)/2 ## as it is either odd-odd or even-even

pjdat[:,left:right] = ff

write_pj_data(pjdat,proj2)

recondets = int(ceil(dts))
reconFOV=detInc*recondets/2.

execute_command("%s %s %s %d %d --offset %f --roi %f,%f,%f,%f" % (PJREC, proj2, variance_iffile, recondets, recondets, corroffset, -reconFOV,reconFOV,-reconFOV,reconFOV))


## 6) fit the real data into the new projections

#pjdat = get_pj_data(proj2)

## 7) reconstruct

## 8) read in reconstructions, output to minc file

## Time to cleanup

# cleanup(workdir)
