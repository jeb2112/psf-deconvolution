/******************************************************************************
 * 
 * libOPT_calcs module -- processing for libOPT python module
 * 
 * Created Aug 6 2006
 * 
 * Copyright 2006, Johnathon R. Walls
 * 
 * ****************************************************************************/


#include <Python.h>
#include <Numeric/arrayobject.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <rfftw.h>

// This line initializes a new python object that will eventually become an exception
static PyObject *libOPTCalcs_Error;
static PyObject *ext_interp2d_by_factor(PyObject *self, PyObject *args);

void initlibOPTcalcs(void);

int main(int argc, char **argv)
{
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initlibOPTcalcs();

    return 1;
}

static PyMethodDef libOPTCalcs_Methods[] = {
	{"interp2d_by_factor", ext_interp2d_by_factor, METH_VARARGS,
		"Performs 2D interpolation on a 2D matrix."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

// Called from main() to initialize everything
void
initlibOPTcalcs(void)
{
    PyObject *m, *d;

	// Initialize the module
    m = Py_InitModule("libOPTcalcs", libOPTCalcs_Methods);
    // Get the dictionary
    d = PyModule_GetDict(m);
    // Create a new type of exception named spam.error
		libOPTCalcs_Error = PyErr_NewException("libOPTcalcs.error", NULL, NULL);
    // Set its value in the dictionary
    PyDict_SetItemString(d, "error", libOPTCalcs_Error);
    import_array();
}
/*
static PyObject *
ext_interp2d(PyObject *self, PyObject *args)
{

	PyArrayObject *mat, *newmat;
	int oldr, oldc, newr, newc, rows, cols;
  float indr,indc,fracr,fracc;
  int findr,findc,i,j;

	if(!PyArg_ParseTuple(args, "iiiiO!O!O!O!", &oldr,&oldc,&newr,&newc,
		&PyArray_Type, &rows, &PyArray_Type, &cols, &PyArray_Type, &mat, &PyArray_Type, &newmat )
		return NULL;

	if(mat->nd != 2) { return NULL; }	
	// Copied from weave code
  for(i=0;i<newr;i++) {
    for(j=0;j<newc;j++) {
      //indr = rows(i);
      indr = *(float *)(rows->data + i*rows->strides[0]);
      //indc = cols(j);
      indc = *(float *)(cols->data + j*cols->strides[0]);
      findr = int(floor(indr));
      findc = int(floor(indc));
      fracr = indr - findr;
      fracc = indc - findc;
      //newmat(i,j) = (1-fracr) * (1-fracc) * mat(findr,findc);
      
      if(findc+1 <= oldc-1)
              newmat(i,j) += (1-fracr) * (fracc) * mat(findr,findc+1);

      if(findr+1 <= oldr-1)
              newmat(i,j) += (fracr) * (1-fracc) * mat(findr+1,findc);

      if((findr+1 <= oldr-1) && (findc+1 <= oldc-1))
              newmat(i,j) +=  (fracr)   * (fracc)   * mat(findr+1,findc+1);
   }
	
}
*/
static PyObject *
ext_interp2d_by_factor(PyObject *self, PyObject *args)
{
	// Honestly, this code blows. Don't use it.  It's so wrong it's not even wrong anymore.
	
	PyArrayObject *mat, *retmat;
	int rowfactor, colfactor, szs[2];
	char *ptr0, *ptr1, *ptr2, *ptr3; 
	int i,j; 
	double rowfrac, colfrac, rowstep, colstep;
	int rowind, colind;

	if(!PyArg_ParseTuple(args, "O!(ii)", &PyArray_Type, &mat, &rowfactor, &colfactor) )
		return NULL;

	if(mat->nd != 2) { return NULL; }

	szs[0] = (mat->dimensions[0] -1) * rowfactor + 1;
	szs[1] = (mat->dimensions[1] -1) * colfactor + 1;
	printf("Old dims: %d x %d.  New sizes: %d x %d\n", mat->dimensions[0], mat->dimensions[1], szs[0],szs[1]);
	
	/* Initialize storage -- if we don't have the memory, bail out now */
	retmat = (PyArrayObject *)PyArray_FromDims(2, szs, PyArray_FLOAT); // returndata->descr->type_num);
	if(retmat == NULL) { return NULL; }	

	rowstep = 1./rowfactor;
	colstep = 1./colfactor;

	for(i=0;i<mat->dimensions[0];i++) {
		for(j=0;j<mat->dimensions[1];j++) {
			printf("At %d,%d we have value %f\n", i, j, *(float *)(mat->data + i*(mat->strides[0]) + j*(mat->strides[1])));
		}
	}

	ptr0 = mat->data;
	ptr2 = retmat->data;
	for(i=0;i<szs[0];i++,ptr2+=retmat->strides[0]) {
		rowind = (int)floor(rowstep*i);
		rowfrac = rowstep*i - rowind;
		ptr0 = mat->data + rowind*(mat->strides[0]);
		ptr3 = ptr2;
		for(j=0;j<szs[1];j++,ptr3+=retmat->strides[1]) {
			colind = (int)floor(colstep*j);
			colfrac = colstep*j - colind;
			ptr1 = ptr0 +  colind*(mat->strides[1]);
			printf("i: %d, j: %d, Rowind: %d Colind: %d\n", i,j, rowind, colind);
			printf("Rowfrac: %f, Colfrac: %f\n", rowfrac, colfrac);
			printf("Ptr1 is %f while ptr3 is %f\n", *(float *)ptr1, *(float *)ptr3);
			*(float *)ptr3 = (1-rowfrac)*(1-colfrac) * *(float *)(ptr1) +
											(1-rowfrac)*( colfrac ) * *(float *)(ptr1+mat->strides[1]) +
											( rowfrac )*(1-colfrac) * *(float *)(ptr1+mat->strides[0]) +
											( rowfrac )*( colfrac ) * *(float *)(ptr1+mat->strides[0]+mat->strides[1]); 
			/*
			(float *)ptr1 = ( (float)((1-rowfactor) * ( 1-colfactor ) *  * (float *)(ptr2 +   (int)rowind   * (mat->strides[0]) +   (int)colind   * (mat->strides[1]))) +
							(float)((1-rowfactor) * (  colfactor  ) *  * (float *)(ptr2 +   (int)rowind   * (mat->strides[0]) + (int)(colind+1) * (mat->strides[1]))) +
							(float)(( rowfactor ) * ( 1-colfactor ) *  * (float *)(ptr2 + (int)(rowind+1) * (mat->strides[0]) +   (int)colind   * (mat->strides[1]))) +
							(float)(( rowfactor ) * (  colfactor  ) *  * (float *)(ptr2 + (int)(rowind+1) * (mat->strides[0]) + (int)(colind+1) * (mat->strides[1]))));
							*/
		}
	}
	
	return retmat;
	
}
