import libOPTcalcs
from scipy import *
b = reshape(arange(16),(4,4)).astype(float32)
c = libOPTcalcs.interp2d_by_factor(b,(1,2))
