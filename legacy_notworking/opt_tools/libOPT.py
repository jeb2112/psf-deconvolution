#!/usr/bin/env python

###################################
#
# libopt.py
# python module to handle OPT tasks
#
# Johnathon R. Walls 2006
#
# TODO:
# - lots
###################################

from scipy import *
import matplotlib
matplotlib.use("Agg")
from pylab import *
import py_minc
import tempfile
import os
import shutil
import commands
from scipy.weave import converters
from scipy.optimize import leastsq
import time
from libOPTconfig import *

## configure libOPT's info prefix
prog_name = sys.argv[0].split("/")[-1]
info_prefix = "### %s informing ###" % prog_name
exec_prefix = "*** %s executing *** " % prog_name

valid_dimensions = ("zspace","yspace","xspace")
valid_slice_calcs = ("sum","std","max")
valid_mincmath_methods = ("add","sub","mult","div","max","maximum","minimum")
valid_minc_datatypes = ("byte","short","int","long","float","double")
valid_recon_filters = ( "abs_bandlimit","abs_sinc","abs_cosine","abs_hamming","abs_hanning","shepp",
    "bandlimit", "sinc", "cosine", "triangle", "hamming", "hanning" )
valid_export_types = ( "pnm" )

## These two should be set according to the os, amd64 or ia32
#if(os.uname()[-1] == "x86_64"):
#    PHM2PJ = "//home/jwalls/Australiawork/packages/CTSim/amd64/tools/ctsimtext phm2pj"
#    PJREC =  "/home/jwalls/Australiawork/packages/CTSim/amd64/tools/ctsimtext pjrec"
#elif(os.uname()[-1] == "i686" or os.uname()[-1] == "i386"):
#    PHM2PJ = "/micehome/jwalls/workspace/CTSim/tools/ctsimtext phm2pj"
#    PJREC =  "/micehome/jwalls/workspace/CTSim/tools/ctsimtext pjrec"

microscope_info = { 0.80 : 0.0125,
                    1.00 : 0.0155,
                    1.30 : 0.0180,
                    1.60 : 0.0225,
                    2.00 : 0.0270,
                    2.50 : 0.0320,
                    3.20 : 0.0390,
                    4.00 : 0.0465,
                    5.00 : 0.0505,
                    6.30 : 0.0620,
                    8.00 : 0.0625,
                    10.0 : 0.0625  }

filter_info = { 'GFP1' : 525,
                'GFP2' : 560,
                'Cy3'  : 590,
                'TXR'  : 620,
                'Rho'  : 620,
                'Cy5'  : 670 }

verbose = 1
debug = 0

if(__name__ == "__main__"):
    print "Don't run me"

## First, some shell utilities

def execute_command(str):
    '''Executes a shell command.  
    Echoes command to stdout if libOPT.verbose is True'''
    if(verbose):
        print "***libOPT Executing***: ", str
    (status, s) = commands.getstatusoutput(str)
    if(status):
        raise SystemError, "Non-zero status output %d and output %s" % (status, s)
    return (status, s)

def provide_feedback(str):
  '''Prints out str to stdout if verbose is set'''
  if(verbose):
    print "###libOPT info###: ", str

## Some file system utilities

def tempfile_handler(dir):
    '''Handles creation of temporary directories or name-generation of temporary files.'''
    try:
        tempfile.tempdir = WORKINGDIR
    except NameError:
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
    '''Generates a temporary file name.'''
    return tempfile_handler(0)

def make_tempdir():
    '''Generates a temporary directory name and creates the directory.'''
    return tempfile_handler(1)

#def cleanup_dir(f):
#    shutil.rmtree(f)
#
#def cleanup_file(f):
#    os.remove(f)

def cleanup(d):
    '''Deletes the file or directory passed as argument.'''
    if(os.path.isdir(d)):
        shutil.rmtree(d)
    else:
        os.remove(d)

## General utilities

def IsEven(d):
    '''Returns 1 if d is even, 0 if it is odd.
    Why does this function exist?'''
    return (not (d%2))

## Numeric utilities

def gaussian_err_calc(parameters, y, x):
    '''Takes the tuple of parameters (a,b,c,d) and returns the error between 
    the data points y and the calculated fit for data points x.'''
    a, b, c, d = parameters
    err = y -  float(c)* exp(-((x-float(a))/float(b))**2) - float(d)
    return err

def exponential_err_calc(parameters, y,x):
    '''Takes the tuple of parameters (a,b,c,d) and returns the error between 
    the data points y and the calculated fit for data points x.'''
    a,b = parameters
    err = y - a*exp(b*x)
    return err

def exponential_fit(vals,ind=False):
    '''Fits the data vals to an exponential curve, returns the parameters as a
    tuple.'''
    if(ind == False):
        ind = arange(len(vals),typecode=float32)
    # initial_guess = [vals[0], max(ind)*log( (vals[0])/(vals[-1]) )]
    initial_guess = [vals[0], ((vals[0]-vals[-1])/vals[0])/len(vals) ]
    if(verbose):
        print "Fitting exponential with an initial guess of", initial_guess
    # print initial_guess, vals, ind
    plsq = leastsq(exponential_err_calc, initial_guess, args=(vals, ind),maxfev=2000)
    
    return plsq[0]
    

def weave_collapse_matrix(retmat, mat, coll):
    code='''
    int i,j,k, c,d,e;
    double val;
    for(i=0;i<nz;i++) {
        for (j=0;j<ny;j++) {
            for(k=0;k<nx;k++) {
                val = 0;
                int loops = 0;
                for(c=0;c<coll;c++) {
                    for(d=0;d<coll;d++) {
                        for(e=0;e<coll;e++) {
                            val += mat(i*coll+c, j*coll+d, k*coll+e);
                            loops++;
                        }
                    }
                }
                if(loops != pow(coll,3)) std::cout << "Mommy!  I did " << loops << " loops!" << std::endl;
                retmat(i,j,k) = val;
            }
        }
    }
    
    '''
    # if(coll == 1): return matr
    (nx,ny,nz) = shape(retmat)
    weave.inline(code,['retmat','mat','coll','nx','ny','nz'],compiler='gcc',type_converters=converters.blitz)
    
def collapse_matrix(matr,coll):
    ''' "Integrates" over the elements in matrix in a collxcollxcoll manner and returns the result. 
        
    So a 20x20 matrix becomes a 5x5 matrix.
    If the larger matrix is not a multiple of the subsampling value, the remainder is ignored.
    Right now, only does a 2 or 3-dim matrix.'''
    # if(coll == 1): return matr
    if(len(shape(matr)) == 3 ):
        shp = (int(floor(shape(matr)[0]/coll)), int(floor(shape(matr)[1]/coll)), int(floor(shape(matr)[2]/coll)))
        retdata = zeros((shp[0],shape(matr)[1],shape(matr)[2]),matr.typecode())
        for i in range(shp[0]):
            retdata[i,:,:] = sum(matr[i*coll:i*coll+coll,:,:],0).astype(matr.typecode()) 
        retdata2 = zeros((shp[0],shp[1],shape(matr)[2]),matr.typecode())
        for j in range(shp[1]):
            retdata2[:,j,:] = sum(retdata[:,j*coll:j*coll+coll,:],1).astype(matr.typecode())
        retdata3 = zeros(shp,matr.typecode())
        del retdata
        for k in range(shp[2]):
            retdata3[:,:,k] = sum(retdata2[:,:,k*coll:k*coll+coll],2).astype(matr.typecode())
    elif(len(shape(matr)) == 2):
        shp = (int(floor(shape(matr)[0]/coll)), int(floor(shape(matr)[1]/coll)))
        retdata2 = zeros((shp[0],shape(matr)[1]),matr.typecode())
        for i in range(shp[0]):
            retdata2[i,:] = sum(matr[i*coll:i*coll+coll,:],0).astype(matr.typecode()) 
        retdata3 = zeros((shp[0],shp[1]),matr.typecode())
        for j in range(shp[1]):
            retdata3[:,j] = sum(retdata2[:,j*coll:j*coll+coll],1).astype(matr.typecode())
    
    return retdata3

def moving_average(a, b):
    '''Averages the input data a with a rect function of pixel length b.'''
    c = ones((b,),float32)
    ret = convolve(a,c,mode=1)
    return ret

def calculate_extents(vals):
    '''Calculates the extent of non-zero information in the array val.
    Returns the beginning pixel and ending pixel.'''
    return (0,len(vals))

def mincmath_handler_extreme(infile1, arg1, outfile, datatype="float", method="mult", dimension="zspace"):
    '''A wrapper function for generic mincmath commands.  The arguments are:
         infile1: must be an existing minc file
         arg1: can be one of const, array, multiple-sliced minc file, single-sliced minc file
         outfile: final file to write out
         datatype: any of the mincmath data types
         method: any of the mincmath methods
         dimension: the dimension to operate along, if arg1 is an array or a single-sliced minc file'''
    
    if(dimension != "zspace" and dimension != "yspace" and dimension != "xspace"):
        raise TypeError, "Dimension %s is not valid" % dimension
    
    f1 = os.path.realpath(infile1)
    (nz, ny, nx) = getminc_dimlengths(f1)
    if(type(arg1) == str): ## It's a file we're talking about
        a1 = os.path.realpath(infile2)
        (nz2, ny2, nx2) = getminc_dimlengths(a1)
    elif(type(arg1) == ArrayType):
        a1 = arg1
    else:
        a1 = float(arg1)
        
    fout = os.path.realpath(outfile)
    
    # First we take care of the case of a 3D minc file multiplied by a constant or a same-sized minc file
    if(type(a1) == float):
        execute_command("mincmath -mult %s -const %s %s" % (f1, str(a1), fout))
        return
    if(type(a1) == str and nz == nz2 and ny == ny2 and nx == nx2):
        execute_command("mincmath -mult %s %s %s" % (f1, str(a1), fout))
        return
    
    ## Next we take care of the case of a 3D minc file multiplied by an array or by a 2D minc file
    
    ## First check that everything is set up properly
    if(type(a1) == ArrayType):
        nelems = len(a1)
        if(dimension == "zspace" and nelems != nz):
            raise SystemError, "The minc file and array do not have the same length along dimension=zspace"
        if(dimension == "yspace" and nelems != ny):
            raise SystemError, "The minc file and array do not have the same length along dimension=yspace"
        if(dimension == "xspace" and nelems != nx):
            raise SystemError, "The minc file and array do not have the same length along dimension=xspace"
    else: ## must be a 2D minc file
        if(dimension == "zspace" and (ny != ny2 or nx != nx2)):
            raise SystemError, "The minc files do not have the same yspace/xspace dimensions (working along zspace)"
        if(dimension == "yspace" and (nz != nz2 or nx != nx2)):
            raise SystemError, "The minc files do not have the same zspace/xspace dimensions (working along yspace)"
        if(dimension == "xspace" and (nz != nz2 or ny != ny2)):
            raise SystemError, "The minc files do not have the same zspace/yspace dimensions (working along xspace)"
    
    ## Ok, we shouldn't need to do more checking.  just do it!
    td = make_tempdir()
    
    if(dimension == "zspace"):
        for i in range(nz):
            if(type(a1) == ArrayType):
                targ = "-const %f" % a1[i]
            else:
                targ = a1
            execute_command("mincreshape -quiet -start %d,0,0 -count 1,%d,%d %s %s/t%04d.mnc" % (i,ny,nx,f1,td,i))
            execute_command("mincmath -quiet -nocheck_dimensions -%s -%s %s/t%04d.mnc %s %s/n%04d.mnc" % (datatype,method,td,i,str(targ),td,i))
            if(verbose):
                print "Completed %d of %d" % (i,nz)
    elif(dimension == "yspace"):
        for i in range(ny):
            if(type(a1) == ArrayType):
                targ = "-const %f" % a1[i]
            else:
                targ = a1
            execute_command("mincreshape -quiet -start 0,%d,0 -count %d,1,%d %s %s/t%04d.mnc" % (i,nz,nx,f1,td,i))
            execute_command("mincmath -quiet -nocheck_dimensions -%s -%s %s/t%04d.mnc %s %s/n%04d.mnc" % (datatype,method,td,i,str(targ),td,i))
            if(verbose):
                print "Completed %d of %d" % (i,ny)
    elif(dimension == "xspace"):
        for i in range(nx):
            if(type(a1) == ArrayType):
                targ = "-const %f" % a1[i]
            else:
                targ = a1
            execute_command("mincreshape -quiet -start 0,0,%d -count %d,%d,1 %s %s/t%04d.mnc" % (i,nz,ny,f1,td,i))
            execute_command("mincmath -quiet -nocheck_dimensions -%s -%s %s/t%04d.mnc %s %s/n%04d.mnc" % (datatype,method,td,i,str(targ),td,i))
            if(verbose):
                print "Completed %d of %d" % (i,nx)

    execute_command("mincconcat -2 -clobber -sequential -%s -concat_dimension %s %s/n*.mnc %s" % (datatype, dimension, td, fout))
    
    cleanup(td)
    
    return

def mincmath_handler(infile1, infile2, outfile, datatype="float", method="mult", dimension="zspace"):
    '''Runs the two minc files through mincmath.  The first file should have multiple slices, the second only one.
      The slices from the first file are extracted using mincreshape, and the calculations performed on the 
      respective slices.  The methods can be any of the mincmath methods, and similarly with datatype.  
      The procedure is done along the dimension dim. sum|std along dimension dim. Default: sum along zspace'''

    f1 = os.path.realpath(infile1)
    f2 = os.path.realpath(infile2)
    fout = os.path.realpath(outfile)
    
    (nz, ny, nx) = getminc_dimlengths(f1)
    
    td = make_tempdir()
    
    if(dimension == "zspace"):
        for i in range(nz):
            execute_command("mincreshape -start %d,0,0 -count 1,%d,%d %s %s/t%04d.mnc" % (i,ny,nx,f1,td,i))
            execute_command("mincmath -nocheck_dimensions -%s -%s %s/t%04d.mnc %s %s/n%04d.mnc" % (datatype,method,td,i,f2,td,i))
    elif(dimension == "yspace"):
        for i in range(ny):
            execute_command("mincreshape -start 0,%d,0 -count %d,1,%d %s %s/t%04d.mnc" % (i,nz,nx,f1,td,i))
            execute_command("mincmath -nocheck_dimensions -%s -%s %s/t%04d.mnc %s %s/n%04d.mnc" % (datatype,method,td,i,f2,td,i))
    else:
        for i in range(nx):
            execute_command("mincreshape -start 0,0,%d -count %d,%d,1 %s %s/t%04d.mnc" % (i,nz,ny,f1,td,i))
            execute_command("mincmath -nocheck_dimensions -%s -%s %s/t%04d.mnc %s %s/n%04d.mnc" % (datatype,method,td,i,f2,td,i))

    execute_command("mincconcat -clobber -2 -sequential -%s -concat_dimension %s %s/n*.mnc %s" % (datatype, dimension, td, fout))
    
    cleanup(td)
    
    return

def calculate_from_minc_slicewise(infile, method="sum", dimension="zspace"):
    '''Uses minc tools to calculate the slicewise sum|std along dimension dim.  Default: sum along zspace'''
    if(not os.path.exists(infile)):
       raise SystemError, "No such file"
    
    (nz, ny, nx) = getminc_dimlengths(infile)
    td = make_tempdir()
    
    if(dimension == "zspace"):
        maxnum = nz
        cmd = "mincreshape -quiet -clobber -start SLICENUM,0,0 -count 1,%d,%d %s %s/slice.mnc" % (ny,nx,infile,td)
    elif(dimension == "yspace"):
        maxnum = ny
        cmd = "mincreshape -quiet -clobber -start 0,SLICENUM,0 -count %d,1,%d %s %s/slice.mnc" % (nz,nx,infile,td)
    elif(dimension == "xspace"):
        maxnum = nx
        cmd = "mincreshape -quiet -clobber -start 0,0,SLICENUM -count %d,%d,1, %s %s/slice.mnc" % (nz,ny,infile,td)
    else:
        raise SystemError, "No such dimension", dimension
    
    vals = zeros((maxnum,),float32)
    for i in range(maxnum):
        ccmd = str(i).join(cmd.split("SLICENUM"))
        execute_command(ccmd)
        ccmd = "mincstats -%s %s/slice.mnc" % (method, td)
        (status, s) = execute_command(ccmd)
        vals[i] = float(s.split('\n')[-1].split(' ')[-1])
        
    cleanup(td)
    
    return vals
    
def calculate_from_minc(infile, method="sum", dimension="zspace", cache="all"):
    '''Uses python arrays to calculate the slicewise sum|std along dimension dim. Default: sum along zspace'''
    mem = py_minc.get_n_bytes_cache_threshold()
    (nz, ny, nx) = getminc_dimlengths(infile)
    if(cache == "slice"):
        print "Slice by slice caching"
        if(dimension == "zspace"):
            py_minc.set_n_bytes_cache_threshold(1*ny*nx*4)
        elif(dimension == "yspace"):
            py_minc.set_n_bytes_cache_threshold(nz*1*nx*4)
        else:
            py_minc.set_n_bytes_cache_threshold(nz*ny*1*4)

    a = openminchandle(infile)
    vals = False

    if(dimension == "zspace"):
        vals = zeros((nz,),float32)
        for i in range(nz):
            print "zspace: %d of %d" % (i,nz)
            if(method == "sum"):
                vals[i] = sum(ravel(a.get_hyperslab((i,0,0),(1,ny,nx))))
            elif(method == "std"):
                vals[i] = std(ravel(a.get_hyperslab((i,0,0),(1,ny,nx))))
            elif(method == "max"):
                vals[i] = max(ravel(a.get_hyperslab((i,0,0),(1,ny,nx))))
    elif(dimension == "yspace"):
        vals = zeros((ny,),float32)
        for i in range(ny):
            if(method == "sum"):
                vals[i] = sum(ravel(a.get_hyperslab((0,i,0),(nz,1,nx))))
            elif(method == "std"):
                vals[i] = std(ravel(a.get_hyperslab((0,i,0),(nz,1,nx))))
            elif(method == "max"):
                vals[i] = max(ravel(a.get_hyperslab((0,i,0),(nz,1,nx))))
    else:
        vals = zeros((nx,),float32)
        for i in range(nx):
            if(method == "sum"):
                vals[i] = sum(ravel(a.get_hyperslab((0,0,i),(nz,ny,1))))
            elif(method == "std"):
                vals[i] = std(ravel(a.get_hyperslab((0,0,i),(nz,ny,1))))
            elif(method == "max"):
                vals[i] = max(ravel(a.get_hyperslab((0,0,i),(nz,ny,1))))
                
    py_minc.set_n_bytes_cache_threshold(mem)
    
    return vals

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
    code=""" printf('hello\n'); """
    weave.inline(code,['oldr','oldc','newr','newc','rows','cols','mat','newmat'],\
                 compiler='gcc',type_converters=converters.blitz)



def interp2d_notworking(mat,sz):
    
    ## Note size can be a multiplier, or a tuple of row/col

    if(len(shape(mat)) != 2):
            print "I can only do 2d matrices"
            return mat

    if(sz == 1):
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

## Normalized 2D cross-correlation routines
## This is done because the 2d cross corr from scipy takes FOREVER

def viewalias(mat, fact):
    '''Interpolates along the first dimension by factor sz.  Does this by zero padding the Fourier space.'''
    if(len(shape(mat)) != 2):
            print "I can only do 2d matrices"
            return mat

    if(fact == 1):
        return mat

    (views, dets) = shape(mat)
    newviews = fact*views
    b = zeros((newviews,dets),float32)
    for i in range(views):
        for j in range(fact):
            ## This is cyclical, so we need to extrapolate from the first view when we hit the last view
            if(i<views-1):
                b[i*fact+j,:] = (1.0*(fact-j)/fact*mat[i,:]+1.0*(j)/fact*mat[i+1,:]).astype(float32)
            else:
                b[i*fact+j,:] = (1.0*(fact-j)/fact*mat[i,:]+1.0*(j)/fact*mat[0,:]).astype(float32)
    
    return b



#    b = fftshift(fft2(fftshift(mat)))
#    c = zeros((shape(b)[0]*fact,shape(b)[1]),Complex64)
#    start = shape(b)[0]/2*fact-shape(b)[0]/2
#    end = start + shape(b)[0]
#    # print start, end, shape(c), shape (c[start:end,:]), shape(b)
#    c[start:end,:] = b
#    d = abs(fftshift(ifft2(fftshift(c))))
#    return d

def SafeDivide(numerator,denominator):
    '''A function which will do the divide without returning NaN.'''
    denom = where(equal(denominator,0),1,denominator)
    numer = where(equal(denominator,0),0,numerator)
    
    return numer / denom

def SafeSqrt(nums):
    '''A function which will do the sqrt without returning imaginary numbers.'''

    newnums = where(less(nums,0),0,nums)

    return sqrt(newnums)

## These functions are duplicated from MatLab normxcorr2

def local_sum(A,m,n):
    '''A sum of some kind.'''
    mA = shape(A)[0]
    nA = shape(A)[1]
    
    B = zeros((2*m + mA - 1, nA * 2 + n -1), float32)
    B[mA:mA+m,n:n+nA] = A.astype(float32);
    s = cumsum(B,0)
    c = s[m:,:] - s[:-m,:]
    s = cumsum(c,1)
    local_sum_A = s[:,n:] - s[:,:-n]
    
    return local_sum_A.astype(float32)

def xcorr2_fast(T,A):
    '''A quick cross correlation in fourier space.'''
    T_size = shape(T)
    A_size = shape(A)
    outsize = zeros(2)
    outsize[0] = A_size[0] + T_size[0] - 1
    outsize[1] = A_size[1] + T_size[1] - 1
    
    aa = freqxcorr2(T,A,outsize);
    
    return aa

def freqxcorr2(a,b,outsize):
    '''The frequency space 2d cross corr'''
    Fa = fft2(rot90(a,2), outsize)
    Fb = fft2(b,outsize)
    xcorr_ab = ifft2(Fa * Fb).real
    return xcorr_ab


def normxcorr2(T,A):
    '''The normalized 2d cross corr copied from Matlab normxcorr2.'''
    xcorr_TA = xcorr2_fast(T,A);
    m = shape(T)[0]
    n = shape(T)[1]
    mn = m*n
    T = T.astype(float32)
    A = A.astype(float32)
    local_sum_A = local_sum(A,m,n)
    local_sum_A2 = local_sum(A*A,m,n)
    ## The following line sometimes returns a ValueError exception
    ## because the value inside sqrt() is negative.
    ## How does this happen?
    tmp = (local_sum_A2 - (local_sum_A ** 2)/mn ) / (mn - 1.)
    denom_A = sqrt(tmp)
    denom_T = std(ravel(T))
    denom = denom_T*denom_A

    numerator = (xcorr_TA - local_sum_A * sum(ravel(T))/mn) / (mn -1)

    C = zeros(shape(numerator))
    ## In the matlab code, only the elements of the denominator
    ## that are non-zero are used to calculate the matrix C
    ## Any values calculated from a 0 denominator value are
    ## set to zero

    ## All zero elements in the denominator become 1
    #denom_new = where(equal(denom,0),1,denom)
    ## All elements in the numerator that would be divided by zero
    ## become zero
    #numerator_new = where(equal(denom,0),0,numerator)

    C = SafeDivide(numerator,denom)


    return C

## Some OPT specific utilities

def get_bg_filename(d):
  '''Parses the dictionary to compile the appropriate bg name'''
  return "%sms_%sx_%s.tif" % (d["Exposure"], d["Magnification"], d["FluorescentFilter"])

def parse_date(s):
  ''' Parses the OPT-prototype date into something more useable'''
  ## Sample format: Wednesday, June 11, 2008
  ## we want returned: june11
  fields = s.split(",")
  if(len(fields) != 3):
    return False
  return ("".join(fields[1].split(" "))).lower()

def parse_annotations(file):
  '''Takes in a ANNOTATIONS_A.txt file from the prototype and returns a dictionary of the values.
  Returns False on a failed read.'''
  ## A sample ANNOTATIONS_A.txt file looks like:
  #@ScanRecord{,
  #% Created by OPT Scanner Script version 1.0
  #      Version = {1.1},
  #% Hardware
  #      MachineID = {HGU01},
  #      Camera = {Retiga EX},
  #      CameraMount = {C-mount 1.0x},
  #      Microscope = {Leica MZFlIII},
  #      Objective = {0.5},
  #      FluorescentSource = {100W mercury lamp},
  #      BrightfieldSource = {20W standard lamp},
  #% Calibration
  #      CalibrationType = {X},
  #      CalibrationRotation = {X.X},
  #      CalibrationAdjust = {X.X},
  #% Experiment
  #      ExperimentID = {538},
  #      Species = {mouse},
  #      Magnification = {3.2},
  #% Specimen
  #      SpecimenID = {dr2},
  #      Strain = {unspecified},
  #      Genotype = {unspecified},
  #      NegativeControl = {NO},
  #% Scan
  #      ScanID = {684},
  #      ScanDate = {Wednesday, June 11, 2008},
  #      ScanTime = {8:01:22 AM},
  #      Username = {jwalls},
  #      RARotation = {CW},
  #% Channel
  #      ChannelID = {1},
  #      ImagingMode = {fluorescent},
  #      BackgroundImages = {NO},
  #      FluorescentLamp = {ON},
  #      BrightfieldLamp = {OFF},
  #      FluorescentFilter = {GFP1},
  #      InfraredFilter = {black},
  #      Exposure = {550},
  #      NumberOfAngles = {400},
  #}
  if(not os.path.exists(file)):
    return False
  a = open(file)
  lines = a.readlines()
  a.close()
  dict = {}
  for line in lines:
    l = line.rstrip().rstrip(",")
    if(l[0] != " "):
      continue
    l = l.lstrip()
    (key, val) = l.split("=")
    key = key.rstrip()
    val = val.lstrip().lstrip("{").rstrip("}")
    dict[key] = val
  
  return dict

def parse_bioptonics_log(file):
  '''Takes in a .log file from the Bioptonics scanner and returns a dictionary of the values.
  Returns False on a failed read.'''
  ## A sample file is:
#[System]
#Scanner=Skyscan3001
#Instrument S/N=13
#Software=Version 1. 3 (build 5)
#Home Directory=C:\Program Files\Bioptonics 3001
#Source Type=GFP1 425/40nm, 475nmLP
#Camera Pixel Size (um)=6.7
#Camera X/Y Ratio=1.0000
#[Acquisition]
#Filter=GFP1 425/40nm, 475nmLP
#Focal Position=182.8
#Optical Axis (line)= 512
#Object to Source (mm)=1000000
#Data Directory=Z:\brain_bioptonics_data\
#Filename Prefix=im
#Number Of Files=  720
#Number Of Rows= 1024
#Number Of Columns= 1024
#Image Pixel Size (um)=11.50
#Image Format=TIFF
#Depth (bits)=16
#Screen LUT=0
#Exposure (ms)=33
#Rotation Step (deg)=0.500
#Use 360 Rotation=YES
#Flat Field Correction=ON
#Rotation Direction=CC
#Type of Detector Motion=STEP AND SHOOT
#Scanning Trajectory=ROUND
#Study Date and Time=May 16, 2008  13:24:26
#[Reconstruction]
#Reconstruction Program=NRecon
#Program Version=Version: 1.5.0
#Program Home Directory=C:\Program Files\Bioptonics 3001
#Dataset Origin=Skyscan3001
#Dataset Prefix=im
#Dataset Directory=Z:\brain_bioptonics_data\
#Time and Date=May 21, 2008  21:27:03
#First Section=0
#Last Section=1022
#Reconstruction duration per slice (seconds)=6.676442
#Postalignment=6.00
#Section to Section Step=1
#Sections Count=1023
#Result File Type=BMP
#Result File Header Length (bytes)=1134
#Result Image Width (pixels)=1024
#Result Image Height (pixels)=1024
#Pixel Size (um)=11.50364
#Reconstruction Angular Range (deg)=360.00
#Use 180+=0
#Angular Step (deg)=0.5000
#Smoothing=0
#Ring Artifact Correction=20
#Draw Scales=ON
#Object Bigger than FOV=OFF
#Reconstruction from ROI=OFF
#Undersampling factor=1
#Threshold for defect pixel mask (%)=50
#CS Static Rotation (deg)=0.0
#Mininum for CS to Image Conversion=0.0021
#Maximum for CS to Image Conversion=0.1270
#HU Calibration=OFF
#BMP LUT=0
#Cone-beam Angle Horiz.(deg)=0.000675
#Cone-beam Angle Vert.(deg)=0.000675
  if(not os.path.exists(file)):
    return False
  a = open(file)
  lines = a.readlines()
  a.close()
  dict = {}
  for line in lines:
    l = line.rstrip()
    if(l[0] == "["):
      continue
    (key, val) = l.split("=")
    key = key.rstrip().lstrip()
    val = val.rstrip().lstrip()
    dict[key] = val
  
  return dict
  

def get_axial_scaling_factor(mag1, mag2, lambda1, lambda2):
    '''Gets the axial scaling factor for resampling based on the ratio of NAs and wavelengths.
    Returns the axial scaling factor necessary for mag2, if mag1 is the base.
    scaling_factor * sample_1 = sample_2, where
    scaling_factor = (lambda_2 * NA_1)/(lambda_1 * NA_2)'''
    if(not microscope_info.has_key(mag1)):
        raise KeyError, 'The microscope config does not have the mag %s in the database' % str(mag1)
    if(not microscope_info.has_key(mag2)):
        raise KeyError, 'The microscope config does not have the mag %s in the database' % str(mag2)
    
    return microscope_info[NA1]/microscope_info[NA2] * lambda2 / lambda1

def get_longitudinal_scaling_factor(NA1, NA2, lambda1, lambda2):
    '''Gets the long. factor for resampling based on the ratio of NAs and wavelengths.
    Returns the long. factor necessary for mag2, if mag1 is the base.
    scaling_factor * sample_1 = sample_2, where
    scaling_factor = (lambda_2 * NA_1**2)/(lambda_1 * NA_2**2)'''
    if(not microscope_info.has_key(mag1)):
        raise KeyError, 'The microscope config does not have the mag %s in the database' % str(mag1)
    if(not microscope_info.has_key(mag2)):
        raise KeyError, 'The microscope config does not have the mag %s in the database' % str(mag2)
    
    return microscope_info[NA1]**2/microscope_info[NA2]**2 * lambda2 / lambda1

## Next, some tiff utilities

def get_tiff_numstacks(filein):
    '''Uses ImageMagick identify to figure out if it's a 2D or 3D tiff file.'''
    cmd = "identify %s" % (filein)
    (stat, out) = execute_command(cmd)
    nz = len(out.split('\n'))
    return nz

def get_tiff_size(filein):
    '''Uses ImageMagick identify to obtain the image size of a 2-D tiff file'''
    cmd = "identify %s" % (filein)
    (stat, out) = execute_command(cmd)
    g = out.split(' ')[2].split('x')
    ret = (int(g[0]),int(g[1]))
    return ret

def get_tiffstack_size(filein):
    '''Uses ImageMagick identify to obtain the image size of a 3-D tiff file'''
    cmd = "identify %s" % (filein)
    (stat, out) = execute_command(cmd)
    ## First we see if it's a TIFF stack
    nz = len(out.split('\n'))
    ## If there is more than one line, then we have a 3D stack
    g = out.split('\n')[0].split(' ')[2].split('x')
    ret = (int(nz), int(g[0]), int(g[1]))
    return ret

## Conversion of tiff to minc file format

def renumber_file(filein, pref):
    '''Converts the silly xaaa.tif to x0000.tif.  returns the new file name.'''
    fnam = filein[len(pref):]
    suff = fnam[3:]
    if(debug):
        print "pref", pref, "fnam", fnam, "suff", suff
    a_num = ord(fnam[0]) - 96
    b_num = ord(fnam[1]) - 96
    c_num = ord(fnam[2]) - 96
    num = (a_num - 1) * 26 * 26 + (b_num-1) * 26 + c_num - 1
    if(debug):
        print "Num", num
    outfile = "%s%04d%s" % (pref, num, suff)
    if(verbose):
        print "Renaming", filein, "to", outfile
    os.rename(filein, outfile)
    return outfile

def tifftominc(filein, fileout, datatype="short", signed=False, invert=False, flip=False, rang=(0,4095)):
    '''Converts a tiff file into a minc file using ImageMagick convert and rawtominc.
    Overwrites any output without checking for clobber.'''
    im_min = rang[0]
    im_max = rang[1]
    imsize = get_tiff_size(filein)
       
    f = "%s.gray" % get_tempfile()
    rf = "%s.flipped" % get_tempfile()
    if(flip):
      cmd = "convert -rotate -90 %s %s" % (filein, rf)
      execute_command(cmd)
      cols = imsize[1]
      rows = imsize[0]
      next = rf
    else:
      cols = imsize[0]
      rows = imsize[1]
      next = filein
    if(signed):
      signed_str = "-signed"
    else:
      signed_str = "-unsigned"
      
    if(invert):
      invert_str = "-negate"
    else:
      invert_str = ""
      
    cmd = 'convert -type Grayscale -size %sx%s %s %s %s' % \
                          (cols, rows, invert_str, next, f)
    execute_command(cmd)
    if(datatype == "byte"):
      cmd = "cat %s | rawtominc -clobber -%s -unsigned -xstep 1 -ystep 1 " \
            "-zstep 1 -ct -origin 0 0 0 -range %f %f -real_range %f %f %s 1 %s %s" % \
            (f, datatype, im_min, im_max, im_min, im_max, fileout, rows, cols)
    else:
      cmd = "cat %s | dd conv=swab | rawtominc -clobber -%s %s -xstep 1 -ystep 1 " \
            "-zstep 1 -ct -origin 0 0 0 -range %f %f -real_range %f %f %s 1 %s %s" % \
            (f, datatype, signed_str, im_min, im_max, im_min, im_max, fileout, rows, cols)
    execute_command(cmd)
    cleanup(f)
    if(os.path.exists(rf)):
      cleanup(rf)

def manytiffs_to_minc(filesin, fileout, datatype="short", signed=False, invert=False, flip=False, rang=(0,4095)):
  '''Converts many tiff files into a single 3D minc file using ImageMagick convert and rawtominc.'''
  
  if(len(rang) != 2):
    raise TypeError, "Range must be a type of two elements"

  td = make_tempdir()  
  for filein in filesin:
    fout = "%s/%s.mnc" % (td, filein.split("/")[-1])
    tifftominc(filein, fout, datatype=datatype, signed=signed, invert=invert, flip=flip, rang=rang)
  
  if(signed):
    signed_str = "-signed"
  else:
    signed_str = "-unsigned"
    
  execute_command("mincconcat -2 -clobber -valid_range %f %f -%s %s -sequential -concat_dimension zspace %s/*.mnc %s" % 
            (rang[0], rang[1], datatype, signed_str, td, fileout))
  
  if(os.path.exists(td)):
    # cleanup(td)
    pass


def split_tiff(filein, odir, prefix="x", hack_endian=False):
  '''Splits a tiff stack into individual tiff files.  Takes into account rotation,
     ranges and the silly endian bug with tiffsplit.'''
 
  infile = os.path.realpath(filein)
  #execute_command('ln -s %s/%s %s/source.tif' % (pwd, filein, td))

  pwd = os.path.realpath('.')
  tf = get_tempfile()
  if(hack_endian):
    execute_command('tiffcp -L %s %s' % (infile, tf))
  else:
    execute_command("cp %s %s" % (infile, tf)) 
  execute_command('tiffsplit %s %s/%s' % (tf, odir, prefix))

  pwd = os.path.realpath('.')
  os.chdir(odir)
  files = os.listdir(odir)
  for file in files:
    if(file[0:len(prefix)] == prefix):
      nf = renumber_file(file, prefix)

  os.chdir(pwd)
  cleanup(tf)
  return


def tiffstack_to_minc(filein, fileout, datatype="short", signed=False, invert=False, flip=False, rang=(0,4095),hack_endian=False):
    '''Converts a tiff file into a minc file using ImageMagick convert and rawtominc.
    Overwrites any output without checking for clobber.'''
    if(len(rang) != 2):
        raise TypeError, "Range must be a tuple of two elements"
    im_min = 0
    im_max = 4095
    imsize = get_tiffstack_size(filein)
    slices = imsize[0]   
    # cols = imsize[1]
    # rows = imsize[2]
    if(signed):
      signed_str = "-signed"
    else:
      signed_str = "-unsigned"
    td = make_tempdir()
    infile = os.path.realpath(filein)
    outfile = os.path.realpath(fileout)
    #execute_command('ln -s %s/%s %s/source.tif' % (pwd, filein, td))

    pwd = os.path.realpath('.')
    os.chdir(td)
    if(hack_endian):
      tf = get_tempfile()
      execute_command('tiffcp -L %s %s' % (infile, tf))
    else:
      tf = infile
    execute_command('tiffsplit %s x' % tf)
    
    files = os.listdir(td)
    for file in files:
        print file
        if(file[0] == 'x'):
            nf = renumber_file(file, 'x')
            tifftominc(nf,"%s.mnc" % nf[0:-4],datatype=datatype,signed=signed,invert=invert,flip=flip,rang=rang)
    os.chdir(pwd)    
    execute_command("mincconcat -2 -valid_range %f %f -%s %s -sequential -concat_dimension zspace %s/x*.mnc %s" % 
                    (rang[0], rang[1], datatype, signed_str, td, outfile))
    
    cleanup(td)
    if(hack_endian):
      cleanup(tf)
  
def export_minc(infile,outdir, export_format="pnm", export_depth="16", dim="zspace",rng = (0,5000)):
  if(len(rng) != 2):
    raise TypeError, "Range must be a tuple of two elements"
  low = float(rng[0])
  high = float(rng[1])
  (nz,ny,nx) = getminc_dimlengths(infile)
  if(dim=="xspace"):
    for i in range(nx):
      execute_command("mincpik -clobber -scale 1 -slice %d -image_range %d %d -coronal -depth %s %s %s/slice%04d.%s" % 
      (i, low, high, export_depth, infile, outdir, i,export_format))
  elif(dim=="yspace"):
    for i in range(ny) :
      execute_command("mincpik -clobber -scale 1 -slice %d -image_range %d %d -sagittal -depth %s %s %s/slice%04d.%s" % 
      (i, low, high, export_depth, infile, outdir, i,export_format))
  else:
    for i in range(nz):
      execute_command("mincpik -clobber -scale 1 -slice %d -image_range %d %d -axial -depth %s %s %s/slice%04d.%s" % 
      (i, low, high, export_depth, infile, outdir, i,export_format))
  
  return
  
  
## Conversion of views to minc file format
## Uses pixel spacing 

def viewtominc(filein, fileout, start, size):
    '''Converts a view into a minc file using ImageMagick convert and rawtominc.
    Differs from tifftominc by setting pixel size and appropriate dimension.
    Does nothing atm.'''
    
    pass

## Now, some minc utilities

def openmincfile(filename, dn=py_minc.FILE_ORDER_DIMENSION_NAMES, tc=int16):
    '''Opens a minc file as ArrayVolume and returns the handle.'''
    a = py_minc.ArrayVolume(filename,dim_names=dn, typecode=tc)
    return a

def openminchandle(filename, dn=py_minc.FILE_ORDER_DIMENSION_NAMES):
    '''Opens a minc file as Volume and returns the handle.'''
    a = py_minc.Volume(filename,dim_names=dn)
    return a

def newmincfile(dims,ord=(py_minc.MIzspace,py_minc.MIyspace,py_minc.MIxspace),
                ncdt=py_minc.NC_FLOAT,tc=float32):
    '''Creates a new minc file and returns the handle.'''
    a = py_minc.ArrayVolume(dims,ord,nc_data_type=ncdt,typecode=tc)
    return a

def getminc_starts(infile):
    '''Retrieves the dimension starts from a minc file.'''
    memback = py_minc.get_n_bytes_cache_threshold()
    py_minc.set_n_bytes_cache_threshold(0)
    a = openminchandle(infile)
    ret = a.get_starts()
    py_minc.set_n_bytes_cache_threshold(memback)
    return ret
    
def getminc_steps(infile):
    '''Retrieves the dimension steps from a minc file.'''
    memback = py_minc.get_n_bytes_cache_threshold()
    py_minc.set_n_bytes_cache_threshold(0)
    a = openminchandle(infile)
    ret = a.get_separations()
    py_minc.set_n_bytes_cache_threshold(memback)
    return ret


def getminc_dimlengths(infile):
    '''Retrieves the dimension lengths from a minc file.'''
    (stat, out) = execute_command("mincinfo -dimlength zspace -dimlength yspace -dimlength xspace %s" % infile)
    sz = out.split('\n')
    ret = (int(sz[0]),int(sz[1]),int(sz[2]))
#    memback = py_minc.get_n_bytes_cache_threshold()
#    py_minc.set_n_bytes_cache_threshold(0)
#    a = py_minc.Volume(infile)
#    sizes = a.get_sizes()
#    del a
#    py_minc.set_n_bytes_cache_threshold(memback)
    return ret

def getminc_history(infile):
    '''Retrieves the history of a minc file and returns it as a list of strings'''
    (stat, out) = execute_command("mincinfo -attvalue :history %s" % infile)
    return out.split('\n')

def getslice(infile, slicenum=350,dimension="zspace"):
    '''Retrieves an image from the slice number slicenum, and returns
    the squeezed array.  Useful for when you simply want just the first image
    and nothing more.'''
    ## This is currently done by setting cache size to zero, then
    ## reading in the appropriate data.  It could also be done using
    ## mincreshape, then reading in that single file.
    memback = py_minc.get_n_bytes_cache_threshold()
    py_minc.set_n_bytes_cache_threshold(0)
    a = py_minc.Volume(infile)
    (nz,ny,nx) = a.get_sizes()
    if(dimension == "zspace"):
        py_minc.set_n_bytes_cache_threshold(ny*nx)
        im = squeeze(a.get_hyperslab((slicenum,0,0),(1,ny,nx)))
    elif(dimension == "yspace"):
        py_minc.set_n_bytes_cache_threshold(nz*nx)
        im = squeeze(a.get_hyperslab((0,slicenum,0),(nz,1,nx)))
    else:
        py_minc.set_n_bytes_cache_threshold(ny*nz)
        im = squeeze(a.get_hyperslab((0,0,slicenum),(nz,ny,1)))
    py_minc.set_n_bytes_cache_threshold(memback)
    del a
    return im

def get_verbose_minc_option(bool=True):
    '''Returns the properly formatted --quiet argument for minc commands.'''
    if(bool):
        return ""
    else:
        return "-quiet"

## Next, some pj/if utilities (file format from CTSim

def create_phantom(str,file):
    '''Creates a .phm file from the supplied phantom string and saves it in
    the specified output file.  Does not check for overwrite.'''
    f = open(file,'w')
    f.writelines(str)
    f.close()

def get_pj_data(file):
    '''Gets the data in a .pj file and returns it as a Numeric array.
    Requires tomography module.'''
    from tomography import projections
    pj = projections.Projections()
    pj.read(file)
    a = pj.get_data()
    return a

def write_pj_data(data,file):
    '''Takes a Numeric array and outputs it to a .pj file.
    Requires tomography module.'''
    from tomography import projections
    pj = projections.Projections()
    pj.read(file)
    pj.set_data(data)
    pj.write(file)

def get_if_data(file):
    '''Gets the data in a .if file and returns it as a Numeric array.
    Requires tomography module.'''
    from tomography import imagefile
    iffile = imagefile.ImageFile()
    iffile.fileRead(file)
    a = iffile.get_data()
    return a

## Core Utilities: These methods are the core functionality of libOPT.py

def correlate_images(im1, im2, blur=1):
    ''' returns a (col,row) tuple of shift necessary for image alignment'''
    
    (cols,rows) = shape(im1)
    
    ## We need a blurring kernel
    x, y = mgrid[-10:10,-10:10]
    r = sqrt(x**2 + y**2)
    lw = blur
    blurf = 1./(2*pi)*exp(-(r**2)/lw**2)

    provide_feedback("Blurring image 1")
    im1 = signal.convolve2d(im1,blurf,'same');
    provide_feedback("Blurring image 2")
    im2 = signal.convolve2d(im2,blurf,'same');

    provide_feedback("Performing cross-correlation")
    my_result = normxcorr2(im1,im2)
    mymax = max(abs(ravel(my_result)))
    listlist = abs(ravel(my_result)).tolist()
    ind = listlist.index(mymax)
    a,b = divmod(ind,rows*2-1)
    col_offset = (a - shape(im1)[0]+1)
    row_offset = (b - shape(im1)[1]+1)
    
    return (col_offset, row_offset)

def get_coarse_jump(infile):
  '''Calculates a coarse jump along the vertical axis between the first and last images.'''
  provide_feedback("Calculating a coarse guess at the vertical jump.  This may take a while ...")
  (views, slices, dets) = getminc_dimlengths(infile)
  im1 = getslice(infile, 0)
  im2 = getslice(infile, views-1)
  return correlate_images(im1, im2, blur=20)[0]

def shift_matrix(m,sh):
  if(sh >= 1.0):
    return m
  if(len(shape(m)) != 2):
    return m
  n = copy(m)
  (nx,ny) = shape(n)
  for i in range(nx-1):
    n[i,:] = ((1-sh)*m[i,:] + sh*m[i+1,:]).astype(float32)
  return n

def correct_jump(infile, outfile, jump=0.0):
  if(jump == 0.0):
    if(verbose):
      provide_feedback("No jump necessary, so the file is unchanged!")
      sys.exit()
  (views, slices, dets) = getminc_dimlengths(infile)
  jump_per_view = jump / float(views)
  execute_command("cp %s %s" % (infile, outfile))
  a = openminchandle(outfile)

  for i in range(views):
    b = a.get_hyperslab((i,0,0),(1,slices,dets))
    b = shift_matrix(b.astype(float32), jump_per_view)
    # a.set_hyperslab((i,0,0),b)
    provide_feedback("Completed view %d of %d" % (i, views))
  return
  
def get_coarse_offset(infile):
    '''Calculates a coarse offset of the rotational axis from the center pixel using normxcorr2.'''
    provide_feedback("Calculating a coarse guess at the alignment.  This may take a while ...")
    (views, slices, dets) = getminc_dimlengths(infile)
    im1 = getslice(infile,0)
    im2 = fliplr(getslice(infile,views/2))
    return correlate_images(im1, im2, blur=12)[1]/2.0
#    (cols,rows) = shape(im1)
#
#    ## We need a blurring kernel
#    x, y = mgrid[-10:10,-10:10]
#    r = sqrt(x**2 + y**2)
#    lw = 8
#    blurf = 1./(2*pi)*exp(-(r**2)/lw**2)
##    blurf = reshape(
##            array([ 0.0001, 0.0011, 0.0029, 0.0011, 0.0001,
##                    0.0011, 0.0215, 0.0585, 0.0215, 0.0011,
##                    0.0029, 0.0585, 0.1592, 0.0585, 0.0029,
##                    0.0011, 0.0215, 0.0585, 0.0215, 0.0011,
##                    0.0001, 0.0011, 0.0029, 0.0011, 0.0001] ),
##                    (5,5))
#    
#    im1 = signal.convolve2d(im1,blurf,'same');
#    im2 = signal.convolve2d(im2,blurf,'same');
#
#    my_result = normxcorr2(im1,im2)
#    mymax = max(abs(ravel(my_result)))
#    listlist = abs(ravel(my_result)).tolist()
#    ind = listlist.index(mymax)
#    a,b = divmod(ind,rows*2-1)
#    col_offset = (a - shape(im1)[0]+1)
#    row_offset = (b - shape(im1)[1]+1)
#    
#    return row_offset/2.0

def get_fine_offset(infile, slice=350, accuracy=0.2, start_offset=False, end_offset=False, plotfile=False ):
    '''Calculates a fine offset of the rotational axis from the center pixel by
    doing many sample recons of various potential offsets.  Takes the max as
    the correct offset, and returns it.'''
    if( (start_offset is False) or (end_offset is False)):
        raise TypeError, "The start_offset or end_offset was not defined"
    
    accuracy = abs(accuracy)

    tf = get_tempfile()
    tdir = make_tempdir()

    (nz,ny,nx) = getminc_dimlengths(infile)
    if(slice > ny):
        raise TypeError, "The stated slice is not within the bounds of slices: %d" % ny
    
    execute_command("mincreshape -clobber -float -start 0,%d,0 -count %d,1,%d %s %s" %
                    (slice, nz,nx, infile, tf))

    first = start_offset
    last = end_offset
    ranges = arange(first,last,accuracy)
    maxdets = int(ceil(nx+2*max(abs(first),abs(last))))
    
    j = 0
    for i in ranges:
        reconstruct(tf, "%s/offset_%04d.mnc" % (tdir, j), offset=i,reconsize=(maxdets,maxdets))
        j += 1

    execute_command("mincconcat -sequential -2 %s -clobber -concat_dimension zspace %s/offset*.mnc %s/all_offsets.mnc" %
                    (get_verbose_minc_option(verbose), tdir, tdir))

    vals = calculate_from_minc("%s/all_offsets.mnc" % tdir, method="std", dimension="zspace")
   
    
    #gplt.plot(vals)
    ml = vals.tolist()
    ind = ml.index(max(vals))
    
    if(plotfile):
      plot(ranges,vals)
      title('std plot of offset cals for %s' % infile)
      xlabel('offset values (%2.2f selected)' % float(ranges[ind]))
      ylabel('std calculation from recon')
      savefig(plotfile)
      
    cleanup(tf)
    cleanup(tdir)
    return(ranges[ind])

def calculate_offset(infile, slice=350, accuracy=0.1, start_offset=False):
    '''Calculates the offset of the rotational axis from the center pixel
    of the minc data set infile.  Assumes views,slices,dets orientation
    according to zspace,yspace,xspace. Slice is the slice to use for testing.
    Accuracy is given in units of pixels. Start_offset skips the first part
    of the code and supplies a near-to-correct value to be further refined.'''

    if(start_offset == False):
        coarse_offset = get_coarse_offset(infile)
    else:
        coarse_offset=start_offset
    if(accuracy < 1.0):
        final_offset = get_fine_offset(infile, slice, accuracy, start_offset=coarse_offset)
    else:
        final_offset = coarse
        
    return final_offset
    
def rebin_mincfile(infile, outfile, binning):
    '''Reads in a minc file, integrates over the pixels to create a smaller array.'''
    a = py_minc.ArrayVolume(infile)
    (null,y,x) = a.get_sizes()
    collapsed = collapse_matrix(squeeze(a.array),binning)

    b = py_minc.ArrayVolume((1,y/binning, x/binning),
                            (py_minc.MIzspace,py_minc.MIyspace,py_minc.MIxspace),
                            typecode=float32,nc_data_type = py_minc.NC_FLOAT)
    b.array[0,:,:] = collapsed.astype(float32)
    b.set_range(float(min(ravel(b.array))),float(max(ravel(b.array))))
    b.output(outfile)

def correct_rolloff_scaling(infile, outfile, filtfile, rolloff=0.3, offset=0,  reconsize = "", pjFOV=1.0, ROI="", offsetview=0, antialias=1, filt="abs_hanning"):
    '''Corrects the loss of peak intensity due to the rolloff filter.
    
    infile must be ordered as views, 1, detectors.
    outfile will be ordered as 1, reconsizexreconsize.'''
    
    ## First we need to generate a phantom to run through the filter
    ## MOST OF THIS CODE SHOULD BE DUPLICATED IN RECONSTRUCT() function
    
    workdir = make_tempdir()
    phantom = "%s/phantom.phm" % workdir
    proj = "%s/projections.pj" % workdir
    proj2 = "%s/projections_scaled.pj" % workdir
    reconim = "%s/outfile.if" % workdir
    reconim2 = "%s/outfile_scaled.if" % workdir
    
    a = openmincfile(infile,tc=float32)
    b = openmincfile(filtfile,tc=float32)
    views, slices, dets = b.get_sizes()
    newviews = views * antialias
    
    detInc = pjFOV*2/dets
    dts = dets + abs(offset*2) + 10
    dets1 = (dts)*sqrt(2)
    FOV1 = pjFOV*dts/dets
    dets2 = int(ceil(dets1))
    
    if((IsEven(dets) and not IsEven(dets2)) or (not IsEven(dets) and IsEven(dets2))):
        dets2 = dets2 + 1
    
    newFOV = FOV1 * dets2/dets1
    
    ## Create the phantom
    ph = "ellipse 0 0 %f %f 0 0\n" % (newFOV, newFOV)
    
    ## At this point we have to insert the many little circles we want to use
    numcircles = int(10.*1./rolloff)-1
    circlespacing = newFOV/numcircles
    for i in range(numcircles):
        ph = ph + "ellipse 0 %f %f %f 0 1\n" % (circlespacing*i, circlespacing/2.0, circlespacing/2.0)
        #ph = ph + "ellipse 0 %f %f %f 0 1\n" % (-circlespacing*i, circlespacing/2.0, circlespacing/2.0)
        #ph = ph + "ellipse %f 0 %f %f 0 1\n" % (circlespacing*i, circlespacing/2.0, circlespacing/2.0)
        #ph = ph + "ellipse %f 0 %f %f 0 1\n" % (-circlespacing*i, circlespacing/2.0, circlespacing/2.0)
                
    create_phantom(ph,phantom)
    ## Simulate projections
    execute_command("%s %s %d %d --offsetview %d --rotangle -1 --phmfile %s" % \
                    (PHM2PJ, proj, dets2, newviews, offsetview, phantom))
    execute_command("cp %s %s" % (proj, proj2))
    
    
    ## Stuff the right data in the projection file
    pjdat = get_pj_data(proj)
    left = (dets2 - dets)/2  ## Note: These will always be even
    right = (dets2 + dets)/2 ## as it is either odd-odd or even-even
    pjdatfft = fftshift(fft2(fftshift(pjdat[:,left:right]))) * squeeze(b.array)
    pjdatback = abs(fftshift(ifft2(fftshift(pjdatfft))))
    pjdat[:,left:right] = pjdatback.astype(float32)
    write_pj_data(pjdat,proj2)

    if(not reconsize):
        nx = int(ceil(dts))
        ny = int(ceil(dts))
    else:
        nx = int(reconsize[0])
        ny = int(reconsize[1])

    if(not ROI):
#        sx=-detInc*nx/2.
#        ex=detInc*nx/2.
#        sy=-detInc*ny/2.
#        ey=detInc*ny/2.
        sx = -newFOV
        ex = newFOV
        sy = -newFOV
        ey = newFOV
    else:
        sx = float(ROI[0])
        ex = float(ROI[1])
        sy = float(ROI[2])
        ey = float(ROI[3])
    
    newdetInc = detInc * dts/nx
    
    execute_command("%s %s %s %d %d --offset %f --roi %f,%f,%f,%f --filter %s" % 
                   (PJREC, proj, reconim, nx, ny, offset, sx,ex,sy,ey, filt))
    execute_command("%s %s %s %d %d --offset %f --roi %f,%f,%f,%f --filter %s" % 
                   (PJREC, proj2, reconim2, nx, ny, offset, sx,ex,sy,ey, filt))
    
    q1 = get_if_data(reconim)
    q2 = get_if_data(reconim2)
    #imshow(q1)
    #imshow(q2)
    r_vals = zeros((numcircles,),float32)
    ys = zeros((numcircles,),float32)
    pts1 = zeros((numcircles,),float32)
    pts2 = zeros((numcircles,),float32)
    
    # print pjFOV, newFOV, ex-sx, detInc, newdetInc
    
    ## Now we calculate the center of each peak
    for i in range(numcircles):
        pt = circlespacing*i
        r_vals[i] = pt
        x_pt = (0 - sx)/newdetInc
        y_pt = (pt- sy)/newdetInc
        y_pt_i = int(floor(y_pt))
        x_pt_i = int(floor(x_pt))
        y_frac = y_pt - y_pt_i
        x_frac = x_pt - x_pt_i
        # print "pt (%f,%f) i (%f,%f) frac (%f,%f)" % (y_pt, x_pt, y_pt_i, x_pt_i, y_frac, x_frac)
        pts1[i]  = (1-y_frac)*(1-x_frac)*q1[ y_pt_i , x_pt_i ]
        pts1[i] += (1-y_frac)*( x_frac )*q1[ y_pt_i ,x_pt_i+1]
        pts1[i] += ( y_frac )*(1-x_frac)*q1[y_pt_i+1, x_pt_i ]
        pts1[i] += ( y_frac )*( x_frac )*q1[y_pt_i+1,x_pt_i+1]
        pts2[i]  = (1-y_frac)*(1-x_frac)*q2[ y_pt_i , x_pt_i ]
        pts2[i] += (1-y_frac)*( x_frac )*q2[ y_pt_i ,x_pt_i+1]
        pts2[i] += ( y_frac )*(1-x_frac)*q2[y_pt_i+1, x_pt_i ]
        pts2[i] += ( y_frac )*( x_frac )*q2[y_pt_i+1,x_pt_i+1]        
        

    #gplt.plot(q2[:,nx/2])
#    gplt.hold("on")
#    gplt.plot(r_vals,pts1,r_vals,pts2,r_vals,pts1/pts2)
#    gplt.plot(arange(nx/2.)/(nx/2.),q2[nx/2:,nx/2]/max(ravel(q2[nx/2:,nx/2])))

    #imshow(q1)
    #imshow(q2)
    #gplt.plot(r_vals, pts, r_vals,pts2,r_vals,pts/pts2)
    
    ## Now to create the spline
    tck = interpolate.splrep(r_vals,pts1/pts2,s=0.001)
    xnew = arange(0,1,0.001)
    ynew = interpolate.splev(xnew,tck,der=0)
    # gplt.plot(xnew,ynew)
    
    final_scale = zeros((ny,nx),float32)
    for i in range(nx):
        x = sx + newdetInc*i
        for j in range(ny):
            y = sy + newdetInc*j
            r = sqrt(x**2+y**2)
            r_i = r/0.001
            frac = r_i - int(r_i)
            #print i,j, r, r_i, len(ynew), ynew[0], shape(final_scale), shape(ynew)
            if(r_i+1 > len(ynew)):
                #print final_scale[i,j], ynew[len(ynew)-1]
                final_scale[i,j] = ynew[len(ynew)-1]
            else:
                final_scale[i,j] = (1-frac)*ynew[int(r_i)] + frac*ynew[int(r_i)+1]
    
    
    
    ggv = newmincfile((1,ny,nx))
    ggv.array[0,:,:] = final_scale.astype(float32)
    ggv.set_range(float(min(ravel(ggv.array))),float(max(ravel(ggv.array))))
    ggv.output(outfile)
    
    cleanup(workdir)
    
def reconstruct(infile, outfile, offset=0, reconsize = "", pjFOV=1.0, ROI="", offsetview=0, antialias=1,
    filt="abs_hanning"):
    '''Performs parallel-ray FBP on a sinogram infile and outputs it to outfile.
    Does not check for overwrite.  Expects that infile exists.
    
    infile must be ordered as a 3D minc file with the dimensions: views, 1, detetectors
    outfile will be ordered as a 3D minc file with dimensions: 1, reconsizeXreconsize
    
    Defaults
    --------
    offset (detector offset of rotational axis): 0
    reconsize (number of pixels to reconstruct, format XxY): dets + 2*offset
    pjFOV (field of view of projections): 1.0
    ROI (what portion of the image to reconstruct): full
    offsetview (offset number of views from zero degrees): 0
    antialias (factor to antialias views) : 1
    filt (filter to use in recon) : abs_hanning
    '''
    ## First set some defaults

    workdir = make_tempdir()
    phantom = "%s/phantom.phm" % workdir
    proj = "%s/projections.pj" % workdir
    reconim = "%s/outfile.if" % workdir
    
    a = openmincfile(infile,tc=float32)
    views, slices, dets = a.get_sizes()
    
    newviews = views * antialias
    
    detInc = pjFOV*2/dets
    dts = dets + abs(offset*2) + 10
    dets1 = (dts)*sqrt(2)
    FOV1 = pjFOV*dts/dets
    dets2 = int(ceil(dets1))
    
    if((IsEven(dets) and not IsEven(dets2)) or (not IsEven(dets) and IsEven(dets2))):
        dets2 = dets2 + 1
    
    newFOV = FOV1 * dets2/dets1
    
    ## Create the phantom
    ph = "ellipse 0 0 %f %f 0 1" % (newFOV, newFOV)
    create_phantom(ph,phantom)
    ## Simulate projections
    execute_command("%s %s %d %d --offsetview %d --rotangle -1 --phmfile %s" % \
                    (PHM2PJ, proj, dets2, newviews, offsetview, phantom))
    
    ## Stuff the right data in the projection file
    pjdat = get_pj_data(proj)
    
    pjdat = zeros(shape(pjdat),float32)
    
    left = (dets2 - dets)/2  ## Note: These will always be even
    right = (dets2 + dets)/2 ## as it is either odd-odd or even-even
    
    #print shape(pjdat[:,left:right]), shape(squeeze(a.array))
    
    ## Now interpolate out the data
    newdat = viewalias(squeeze(a.array),antialias)
    
    pjdat[:,left:right] = newdat.astype(float32)
    
    write_pj_data(pjdat,proj)
    
    ## Have to do some more initializing
    
    if(not reconsize):
        nx = int(ceil(dts))
        ny = int(ceil(dts))
    else:
        nx = int(reconsize[0])
        ny = int(reconsize[1])
        
    if(not ROI):
        sx=-detInc*nx/2.
        ex=detInc*nx/2.
        sy=-detInc*ny/2.
        ey=detInc*ny/2.
    else:
        sx = float(ROI[0])
        ex = float(ROI[1])
        sy = float(ROI[2])
        ey = float(ROI[3])
    
    ## Reconstruct the projection file
    cmd = "%s %s %s %d %d --offset %f --roi %f,%f,%f,%f --filter %s" %  \
          (PJREC, proj, reconim, nx, ny, offset, sx,ex,sy,ey, filt)
    execute_command(cmd)
    
    ifdata = get_if_data(reconim)
    
    b = newmincfile((slices,ny,nx))
    b.array[0,:,:] = ifdata
    b.set_range(float(min(ravel(ifdata))),float(max(ravel(ifdata))))
    ## Now pull in history and set
    hist = getminc_history(infile)
    if(hist == None):
        hist = [ "%s>>>%s" % (time.ctime(), ' '.join(sys.argv)) ]
    else:
        hist.append("%s>>>%s" % (time.ctime(), ' '.join(sys.argv)))
    b.output(outfile, '\n'.join(hist))
    if(verbose):
        print "**libOPT** : outputted reconstructed file to", outfile
    cleanup(workdir)

