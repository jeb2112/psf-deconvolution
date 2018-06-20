/***********************************
 *
 * libfdp.c - takes files in, 3DFFT, filters them
 * Johnathon Walls 2005
 *
 ***********************************/

#include "fdr.h"

int verbose = 2;
int use_fftshift = 1;

/***************************************************************
 ***************************************************************
 ******************* UTILITY FUNCTIONS *************************
 ***************************************************************
 ***************************************************************/

int print_data(long n3, long n2, long n1, fdr_complex* data) {

  long k,j,i, index;
  for(k=0;k<n3;k++) {
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		fprintf(stdout, "Z=%d Y=%d X=%d %f + %fi\n", k,j,i,data[index].re,data[index].im);
	  }
	}
  }
  return 1;
}


float sgsum(int elems, float *data) 
{
  int i;
  float sm = 0;
  printf("Here goes, %d elems\n",elems);
  for (i=0;i<elems;i++) {
	//if(i%1000==5)
	//  printf("%d %f\n", i,data[i]);
	sm += data[i];
  }
  return sm;
}

float mymax(int elems, float *data) {
  int i;
  float mx = 0;
  for(i=0;i<elems;i++) {
	if(data[i] > mx) 
	  mx = data[i];
  }
  return mx;
}

float mymin(int elems, float *data) {
  int i;
  float mn = 999999999999.;
  for(i=0;i<elems;i++) {
	if(data[i] < mn) 
	  mn = data[i];
  }
  return mn;
}

float absmax(int elems,fdr_complex *data) {
  int i;
  float mx = 0;
  for(i=0;i<elems;i++) 
	if(sqrt(pow(data[i].re,2)+pow(data[i].im,2)) > mx)
	  mx = sqrt(pow(data[i].re,2)+pow(data[i].im,2));
  return mx;
}

float complexmax(int elems,fdr_complex *data) {
  int i;
  float mx = 0;
  for(i=0;i<elems;i++) 
	if(data[i].re > mx)
	  mx = data[i].re;
  return mx;
}

float complexmin(int elems,fdr_complex *data) {
  int i;
  float mx = 999999999999;
  for(i=0;i<elems;i++) 
	if(data[i].re < mx)
	  mx = data[i].re;
  return mx;
}

int correct_fft_phaseshift(long Z, long Y, long X, float *data) {

  long i,j,k, index;

  for(k=0;k<Z;k++) {
	for(j=0;j<Y;j++) {
	  for(i=0;i<Z*Y*X;i++) {
		index = k*Y*X + j*X + i;
		data[index] *= pow((-1),(index));
	  }
	}
  }
  return 1;
}
	
int shift_complex_to_float(long elems, float *data) {

  fdr_complex *tmp;
  long i;

  tmp = (fdr_complex *)data;
  for(i=0;i<elems;i++,tmp++) {
	data[i] = (*tmp).re;
  }

  return 1;
}
	
int shift_float_to_complex(long elems, fdr_complex *data) {
  float *tmp;
  long i;

  tmp = (float *)data;
  tmp += elems-1;
  // Start at the end of the data and transpose
  for(i=elems-1;i>=0;i--) 
  {
		//fprintf(stdout,"Index is %d\n",i);
		data[i]= (float) *tmp + 0*I;
		//fprintf(stdout,"pointer value is %f, data is %f %f\n", *tmp, data[i].re, data[i].im);
		tmp = tmp - 1;
  }

  return 1;
}

int calculate_abs(long nelems, fdr_complex *data) {

  long i;
  fdr_complex *tmp;
  float *tmp2;
  float x;
  tmp = (fdr_complex *)data;
  tmp2 = (float *)data;
  /* Do not parallelize this with pragma omp as it causes syncing issues */
  //#pragma omp parallel for
  for(i=0;i<nelems;i++) {
	x = sqrt( tmp[i].re*tmp[i].re + tmp[i].im*tmp[i].im );
	tmp2[i] = x;
	if(isinf(tmp2[i])) 
	  fprintf(stdout, "Infinity calculated at index %d, from %f %f\n", i, tmp[i].re, tmp[i].im);
	//*tmp2 = sqrt( (*tmp).re * (*tmp).re + (*tmp).im * (*tmp).im );
	//	fprintf(stdout, "Index %03d re %f im %f abs %f\n", i,(*tmp).re, (*tmp).im, tmp2);
	//data[i] = tmp2;
	//tmp = tmp + 1;
	//tmp2 = tmp2 + 1;
  }

  return 1;
}

int calculate_phase(long nelems, fdr_complex *data) {
  long i;
  fdr_complex *tmp;
  float *tmp2;
  float x;
  tmp = (fdr_complex *)data;
  tmp2 = (float *)data;
  /* Do not parallelize this with pragma omp as it causes syncing issues */
  //#pragma omp parallel for
  for(i=0;i<nelems;i++) {
	// x = sqrt( tmp[i].re*tmp[i].re + tmp[i].im*tmp[i].im );
	x = atan2f(tmp[i].im, tmp[i].re);
	tmp2[i] = x;
	if(isinf(tmp2[i])) 
	  fprintf(stdout, "Infinity calculated at index %d, from %f %f\n", i, tmp[i].re, tmp[i].im);
	//*tmp2 = sqrt( (*tmp).re * (*tmp).re + (*tmp).im * (*tmp).im );
	//	fprintf(stdout, "Index %03d re %f im %f abs %f\n", i,(*tmp).re, (*tmp).im, tmp2);
	//data[i] = tmp2;
	//tmp = tmp + 1;
	//tmp2 = tmp2 + 1;
  }

  return 1;
}


/***************************************************************
 ***************************************************************
 ************** MNC FILE MANIPULATION FUNCTIONS ****************
 ***************************************************************
 ***************************************************************/


int open_minc_file(char *filein, long *nz, long *ny, long *nx, void **data, int flags) 
{

  mihandle_t in_volume;
  midimhandle_t in_dimensions[3];
  unsigned int in_dimsizes[3];
  unsigned long start[3], count[3],nelems;

  int result;

  /* Get volume handle from minc file */
  result = miopen_volume(filein,MI2_OPEN_READ, &in_volume);
  if (result == MI_ERROR) {
	fprintf(stderr, "Error opening the input file %s.\n", filein);
	return 0;
  } 
  else if(verbose>1) fprintf(stdout,"Opened file\n");

  /* Get the dimension handles */
  result = miget_volume_dimensions (in_volume, MI_DIMCLASS_SPATIAL, 
									MI_DIMATTR_ALL, MI_DIMORDER_FILE, 
									3, in_dimensions);
  if (result == MI_ERROR) {
    fprintf(stderr, "Error getting dimension handles %d.\n", result);
	return 0;
  }
  else if(verbose>1) fprintf(stdout,"Got dimension handles\n");

  /* Get the dimension sizes */
  result = miget_dimension_sizes (in_dimensions,3,in_dimsizes);
  if (result == MI_ERROR) {
    fprintf(stderr, "Error getting dimension sizes.\n");
	return 0;
  }
  else if(verbose>1) fprintf(stdout,"Got the dimension sizes Z=%d Y=%d X=%d\n", 
					  in_dimsizes[0],in_dimsizes[1],in_dimsizes[2]);
  nelems = in_dimsizes[0]*in_dimsizes[1]*in_dimsizes[2];

  /* Allocate the memory */
  if(flags == REAL_AS_COMPLEX || flags == COMPLEX_AS_COMPLEX || flags == COMPLEX_AS_REAL) 
	*data = (fdr_complex *)calloc(nelems,sizeof(fdr_complex));
  else
	*data = (float *)calloc(nelems,sizeof(float));
  if(*data == NULL) {
	fprintf(stderr, "Could not allocate memory!\n");
	return 0;
  }
  else if(verbose>1) fprintf(stdout,"Allocated memory\n");

  /* Retrieve the data from the file */
  start[0] = start[1] = start[2] = 0;
  count[0] = (unsigned long)in_dimsizes[0];
  count[1] = (unsigned long)in_dimsizes[1];
  count[2] = (unsigned long)in_dimsizes[2];
  if(flags == REAL_AS_COMPLEX || flags == REAL_AS_REAL)
	result = miget_voxel_value_hyperslab(in_volume,MI_TYPE_FLOAT,start,count, (float *) *data);
  else if(flags == COMPLEX_AS_COMPLEX || flags == COMPLEX_AS_REAL) 
	result = miget_voxel_value_hyperslab(in_volume,MI_TYPE_FCOMPLEX,start,count,(fdr_complex *) *data);
  if(result == MI_ERROR) {
	fprintf(stderr, "Error getting data.\n");
	return 0;
  }
  else if(verbose>1) fprintf(stdout,"Retrieved data\n");

  /* Perform any necessary conversions */
  if(flags == REAL_AS_COMPLEX) {
	shift_float_to_complex(nelems, (float *) *data);
	if(verbose>1) fprintf(stdout, "Rearranged data from float to complex\n");
  }
  else if(flags == COMPLEX_AS_REAL) {
	calculate_abs(nelems, (fdr_complex *) *data);
	if(verbose) fprintf(stdout, "Rearranged data from complex to float\n");
  }

  /* Close the file handle, we've got what we came for */
  miclose_volume(in_volume);

  /* Set the variables whose value was requested*/
  *nz = in_dimsizes[0];
  *ny = in_dimsizes[1];
  *nx = in_dimsizes[2];
  if(verbose>1) fprintf(stdout, "Completed open_minc\n");

  return 1;

}

int get_minc_subvol(char *filein, unsigned long *start, unsigned long *count, long *nz, long *ny, long *nx, void **data, int flags) 
{

  mihandle_t in_volume;
  midimhandle_t in_dimensions[3];
  unsigned int in_dimsizes[3];
  unsigned long nelems;

  int result;

  /* Get volume handle from minc file */
  result = miopen_volume(filein,MI2_OPEN_READ, &in_volume);
  if (result == MI_ERROR) {
	fprintf(stderr, "Error opening the input file %s.\n", filein);
	return 0;
  } 
  else if(verbose>1) fprintf(stdout,"Opened file\n");

  /* Get the dimension handles */
  result = miget_volume_dimensions (in_volume, MI_DIMCLASS_SPATIAL, 
									MI_DIMATTR_ALL, MI_DIMORDER_FILE, 
									3, in_dimensions);
  if (result == MI_ERROR) {
    fprintf(stderr, "Error getting dimension handles %d.\n", result);
	return 0;
  }
  else if(verbose>1) fprintf(stdout,"Got dimension handles\n");

  /* Get the dimension sizes */
  result = miget_dimension_sizes (in_dimensions,3,in_dimsizes);
  if (result == MI_ERROR) {
    fprintf(stderr, "Error getting dimension sizes.\n");
	return 0;
  }
  else if(verbose>1) fprintf(stdout,"Got the dimension sizes Z=%d Y=%d X=%d\n", 
					  in_dimsizes[0],in_dimsizes[1],in_dimsizes[2]);

  /* Verify that we're not asking the impossible */
  if(start[0] + count[0] > in_dimsizes[0] || start[1] + count[1] > in_dimsizes[1] || start[2] + count[2] > in_dimsizes[2]) {
    fprintf(stderr, "Insufficient dimension space!  (%d,%d,%d)+(%d,%d,%d) > (%d,%d,%d)\n", start[0],start[1],start[2],
            count[0],count[1],count[2],in_dimsizes[0],in_dimsizes[1],in_dimsizes[2]);
    return 0;
  }
  nelems = count[0]*count[1]*count[2];

  /* Allocate the memory */
  if(flags == REAL_AS_COMPLEX || flags == COMPLEX_AS_COMPLEX || flags == COMPLEX_AS_REAL) 
	*data = (fdr_complex *)calloc(nelems,sizeof(fdr_complex));
  else
	*data = (float *)calloc(nelems,sizeof(float));
  if(*data == NULL) {
	fprintf(stderr, "Could not allocate memory!\n");
	return 0;
  }
  else if(verbose>1) fprintf(stdout,"Allocated memory\n");


  /* Retrieve data from file */
  if(flags == REAL_AS_COMPLEX || flags == REAL_AS_REAL)
	result = miget_voxel_value_hyperslab(in_volume,MI_TYPE_FLOAT,start,count, (float *) *data);
  else if(flags == COMPLEX_AS_COMPLEX || flags == COMPLEX_AS_REAL) 
	result = miget_voxel_value_hyperslab(in_volume,MI_TYPE_FCOMPLEX,start,count,(fdr_complex *) *data);
  if(result == MI_ERROR) {
	fprintf(stderr, "Error getting data.\n");
	return 0;
  }
  else if(verbose>1) fprintf(stdout,"Retrieved data\n");

  /* Perform any necessary conversions */
  if(flags == REAL_AS_COMPLEX) {
	shift_float_to_complex(nelems, (float *) *data);
	if(verbose>1) fprintf(stdout, "Rearranged data from float to complex\n");
  }
  else if(flags == COMPLEX_AS_REAL) {
	calculate_abs(nelems, (fdr_complex *) *data);
	if(verbose) fprintf(stdout, "Rearranged data from complex to float\n");
  }

  /* Close the file handle, we've got what we came for */
  miclose_volume(in_volume);

  /* Set the variables whose value was requested*/
  *nz = in_dimsizes[0];
  *ny = in_dimsizes[1];
  *nx = in_dimsizes[2];
  if(verbose>1) fprintf(stdout, "Completed open_minc\n");

  return 1;

}


int get_minc_dimensions_from_handle(mihandle_t filein, long *nz, long *ny, long *nx)
{
  midimhandle_t in_dimensions[3];
  unsigned int in_dimsizes[3];
  int result;


  /* Get the dimension handles */
  result = miget_volume_dimensions (in_volume, MI_DIMCLASS_SPATIAL, 
									MI_DIMATTR_ALL, MI_DIMORDER_FILE, 
									3, in_dimensions);
  if (result == MI_ERROR) 
  {
    fprintf(stderr, "Error getting dimension handles %d.\n", result);
		return 0;
  }
  else if(verbose>1) fprintf(stdout,"Got dimension handles\n");

  /* Get the dimension sizes */
  result = miget_dimension_sizes (in_dimensions,3,in_dimsizes);
  if (result == MI_ERROR) 
  {
    fprintf(stderr, "Error getting dimension sizes.\n");
		return 0;
  }
  else if(verbose>1) 
  	fprintf(stdout,"Retrieved dimensions Z=%d, Y=%d, X=%d\n",
  									in_dimsizes[0],in_dimsizes[1],in_dimsizes[2]);

  /* Set the variables we requested */
  *nz = in_dimsizes[0];
  *ny = in_dimsizes[1];
  *nx = in_dimsizes[2];

  return 1;
}
	

int get_minc_dimensions(char *filein, long *nz, long *ny, long *nx)
{
  mihandle_t in_volume;
  midimhandle_t in_dimensions[3];
  unsigned int in_dimsizes[3];

  int result;

  /* Get volume handle from minc file */
  result = miopen_volume(filein,MI2_OPEN_READ, &in_volume);
  if (result == MI_ERROR) 
  {
		fprintf(stderr, "Error opening the input file %s.\n", filein);
		return 0;
  } 
  else if(verbose>1) 
  	fprintf(stdout,"Opened file\n");

  /* Get the dimension handles */
  result = get_minc_dimensions_from_handle(in_volume, nz, ny, nx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout,"Got dimension handles\n");

  /* Close the volume, we're done here */
  miclose_volume(in_volume);

  return 1;

}

int write_minc_file(char *fileout, long nz, long ny, long nx, void *data, int flags)
{
  mihandle_t out_volume;
  midimhandle_t out_dimensions[3];
  unsigned int out_dimsizes[3];
  unsigned long out_count[3];
  unsigned long start[3];
  long nelems;

  float tmp, tmp2;
  int result;
  
  if(verbose>1) fprintf(stdout, "Entering write minc file\n");
  
  start[0] = start[1] = start[2] = 0;
  out_dimsizes[0] = nz;
  out_dimsizes[1] = ny;
  out_dimsizes[2] = nx;
  nelems = nz*ny*nx;

  /* Create the three dimensions dimension */
  result = micreate_dimension("xspace",MI_DIMCLASS_SPATIAL,MI_DIMATTR_REGULARLY_SAMPLED,
							  out_dimsizes[2],&out_dimensions[2]);
  if(result == MI_ERROR) { 
	fprintf(stderr, "Error creating output dimensions\n"); 
	return 0;
  }  
  result = micreate_dimension("yspace",MI_DIMCLASS_SPATIAL,MI_DIMATTR_REGULARLY_SAMPLED,
							  out_dimsizes[1],&out_dimensions[1]);
  if(result == MI_ERROR) { 
	fprintf(stderr, "Error creating output dimensions\n"); 
	return 0;
  }  
  result = micreate_dimension("zspace",MI_DIMCLASS_SPATIAL,MI_DIMATTR_REGULARLY_SAMPLED,
							  out_dimsizes[0],&out_dimensions[0]);
  if(result == MI_ERROR) { 
	fprintf(stderr, "Error creating output dimensions\n"); 
	return 0;
  }  
  else if(verbose>1) fprintf(stdout, "Created dimensions for output file.\n");

  /* Create the file on the disk and get volume handle*/
  if(flags == COMPLEX_AS_COMPLEX || flags == REAL_AS_COMPLEX) {
	result = micreate_volume(fileout, 3, out_dimensions,MI_TYPE_FCOMPLEX, MI_CLASS_COMPLEX,NULL,&out_volume);
  }
  else if(flags == COMPLEX_AS_REAL || flags == REAL_AS_REAL)
	result = micreate_volume(fileout, 3, out_dimensions,MI_TYPE_FLOAT, MI_CLASS_REAL,NULL,&out_volume);
  if(result == MI_ERROR) {
	fprintf(stderr, "Error creating output volume\n");
	return 0;
  }
  else if(verbose>1) fprintf(stdout, "Created volume for output.\n");

  /* Perform any conversions */
  if(flags == COMPLEX_AS_REAL) {
	calculate_abs(nelems,(fdr_complex *)data);
	if(verbose) fprintf(stdout, "Calculated magnitude data\n");
	if(verbose) 
	  fprintf(stdout, "****NOTE**** Data has been changed, any operations must take this into account!\n");
  }
  else if(flags == REAL_AS_COMPLEX) {
	shift_float_to_complex(nelems, (float *)data);
	if(verbose) fprintf(stdout, "Rearranged data from float to complex\n");
	if(verbose) 
	  fprintf(stdout, "****NOTE**** Data has been changed, any operations must take this into account!\n");
  }

  /* prep the count due to long/int mismatch */
  out_count[0] = out_dimsizes[0];
  out_count[1] = out_dimsizes[1];
  out_count[2] = out_dimsizes[2];
  if(verbose>1) fprintf(stdout, "Output dimension sizes Z=%d Y=%d X=%d\n", 
					 (int)out_count[0], (int)out_count[1], (int)out_count[2]);

  /* Create the image associated with the volume*/
  result = micreate_volume_image(out_volume);
  if(result == MI_ERROR) {
	fprintf(stderr, "Error creating output image\n");
	return 0;
  }
  else if(verbose>1) fprintf(stdout, "Created output image for volume.\n");

  /* Get the minimum and maximum and set the file's range */
  if(flags == REAL_AS_REAL || flags == COMPLEX_AS_REAL) {
	tmp = mymax(nelems, (float *) data);
	tmp2 = mymin(nelems, (float *) data);
  }
  else if(flags == REAL_AS_COMPLEX || flags == COMPLEX_AS_COMPLEX) {
	tmp = complexmax(nelems, (fdr_complex*) data);
	tmp2 = complexmin(nelems, (fdr_complex*) data);
  }
  result = miset_volume_valid_range(out_volume, tmp, tmp2);
  if(result == MI_ERROR) {
	fprintf(stderr, "Could not set valid range\n");
	return 0;
  }
  if(verbose>1) fprintf(stdout, "Set output valid ranges\n");

  /* If the data is real, also set the volume range - this doesn't work for complex data*/
  if(flags == COMPLEX_AS_REAL || flags == REAL_AS_REAL) {
	result = miset_volume_range(out_volume, tmp, tmp2);
	if(result == MI_ERROR) {
	  fprintf(stderr, "Could not set valid range 2\n");
	  return 0;
	}
	if(verbose>1) fprintf(stdout, "Set output valid ranges 2\n");
  }

  /* Write out the data */
  if(flags == REAL_AS_COMPLEX || flags == COMPLEX_AS_COMPLEX)
	result = miset_voxel_value_hyperslab(out_volume,MI_TYPE_FCOMPLEX,start,out_count,(mifcomplex_t *)data);
  else if(flags == COMPLEX_AS_REAL || flags == REAL_AS_REAL)
	result = miset_voxel_value_hyperslab(out_volume,MI_TYPE_FLOAT, start, out_count, (float *)data);
  if(result == MI_ERROR) {
	fprintf(stderr, "Error setting output data.\n");
	return 0;
  }

  /* Close the volume, we've completed what we want */
  miclose_volume(out_volume);

  if(verbose>1) fprintf(stdout, "Set output hyperslab\n");

  return 1;

}

int cfftshift(long Z,long Y,long X, fdr_complex* data) {

  fdr_complex tmp;
  long i,j,k, index1, index2, Y2, Z2, X2;;
  
  if (data == NULL) {
	fprintf(stderr, "Fftshift called with no data\n");
	return 0;
  }

  Z2 = Z/2;
  Y2 = Y/2;
  X2 = X/2;
#pragma omp parallel for private(index1,index2,tmp,i,j)
  for(k=0;k<Z2; k++) {
	// fprintf(stdout, "FFTSHIFT Z=%d\n", k);
	for(j=0;j<Y2;j++) {
	  for(i=0;i<X2;i++) {

		// Swap quadrant 1 and 5
		index1 = k*Y*X + j*X + i;
		index2 = (k+Z2)*Y*X + (j+Y2)*X + i+X2;
		tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
		// Swap quadrant 2 and 6
		index1 = k*Y*X + j*X + i+X2;
		index2 = (k+Z2)*Y*X + (j+Y2)*X + i;
		tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
		// Swap quadrant 3 and 7
		index1 = k*Y*X + (j+Y2)*X + i;
		index2 = (k+Z2)*Y*X + j*X + i+X2;
		tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
		// Swap quadrant 4 and 8
		index1 = k*Y*X + (j+Y2)*X + i+X2;
		index2 = (k+Z2)*Y*X + j*X + i;
		tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
	  }
	}
  }

  return 1;
}

/***************************************************************
 ***************************************************************
 ********************** FFT FUNCTIONS **************************
 ***************************************************************
 ***************************************************************/



int c2dfftshift(long Y,long X, fdr_complex* data) {

  fdr_complex tmp;
  long i,j,k, index1, index2, Y2, X2;
  
  if (data == NULL) {
	fprintf(stderr, "Fftshift called with no data\n");
	return 0;
  }

  Y2 = Y/2;
  X2 = X/2;
#pragma omp parallel for private(index1,index2,tmp,i,j)  
  for(j=0;j<Y2;j++) {
	for(i=0;i<X2;i++) {
	  
	  // Swap quadrant 1 and 3
	  index1 = j*X + i;
	  index2 = (j+Y2)*X + i+X2;
	  tmp = data[index1];
	  data[index1] = data[index2];
	  data[index2] = tmp;
	  // Swap quadrant 2 and 4
	  index1 = j*X + i+X2;
	  index2 = (j+Y2)*X + i;
	  tmp = data[index1];
	  data[index1] = data[index2];
	  data[index2] = tmp;
	}
  }
  return 1;
}

int fftshift(long Z,long Y,long X, float* data) {

  float tmp;
  long i,j,k, index1, index2, Y2, Z2, X2;;
  
  if (data == NULL) {
	fprintf(stderr, "Fftshift called with no data\n");
	return 0;
  }

  Z2 = Z/2;
  Y2 = Y/2;
  X2 = X/2;
  for(k=0;k<Z2; k++) {
	// fprintf(stdout, "FFTSHIFT Z=%d\n", k);
	for(j=0;j<Y2;j++) {
	  for(i=0;i<X2;i++) {

		// Swap quadrant 1 and 5
		index1 = k*Y*X + j*X + i;
		index2 = (k+Z2)*Y*X + (j+Y2)*X + i+X2;
		tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
		// Swap quadrant 2 and 6
		index1 = k*Y*X + j*X + i+X2;
		index2 = (k+Z2)*Y*X + (j+Y2)*X + i;
		tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
		// Swap quadrant 3 and 7
		index1 = k*Y*X + (j+Y2)*X + i;
		index2 = (k+Z2)*Y*X + j*X + i+X2;
		tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
		// Swap quadrant 4 and 8
		index1 = k*Y*X + (j+Y2)*X + i+X2;
		index2 = (k+Z2)*Y*X + j*X + i;
		tmp = data[index1];
		data[index1] = data[index2];
		data[index2] = tmp;
	  }
	}
  }

  return 1;
}


int sg_ccfft3d(long n3, long n2, long n1, fdr_complex *data, int direction)
{
  int result;
  long ldx, ldx2, ldy, ldy2;
  long maxdim;
  float scale;
  //long i,j,k;

  float *table = NULL;
  float *work = NULL;
  int isys[2];

  if(data == NULL) {
	fprintf(stderr, "sg_ccfft3d called with NULL pointers\n");
	return 0;
  }
  
  /* Set up coordinates */
  ldx = n1;
  ldx2 = n2;
  ldy = n1;
  ldy2 = n2;
  maxdim = n1>n2 ? (n1>n3 ? n1:n3) : (n2>n3 ? n2:n3);
  
  /* Allocate temporary storage memory */
  table = (float *)calloc(((2*n1+256) + (2*n2+256) + (2*n3+256)),sizeof(float));
  work = (float *)calloc((2*maxdim),sizeof(float));

  if(table == NULL || work == NULL) {
	fprintf(stderr, "Could not allocate memory!\n");
	return 0;
  }

  if(verbose>1) fprintf(stdout,"Allocated memory\n");

  if(use_fftshift)
	cfftshift(n3,n2,n1,data);

  if(verbose>1) fprintf(stdout,"Completed fftshift\n");

  /* Perform 3D FFT using SGI's SCSL libraries */
  isys[0] = 1;
  if(direction == FORWARD_FFT) 
	scale = 1.0f;
  else if(direction == INVERSE_FFT)
	scale = 1.0f/(1.0*n1*n2*n3);
  else {
	fprintf(stderr, "FFT DIRECTION IS NOT DEFINED, ASSUMING FORWARD\n");
	scale = 1.0f;
  }
  /* in-place 3dFFTs */
  ccfft3d(0, n1, n2, n3, scale, 
		  (fdr_complex *)data, ldx, ldx2, 
		  (fdr_complex *)data, ldy, ldy2,
		  table, work, isys);
  // if(verbose>1) fprintf(stdout, "Prepped table\n");
  ccfft3d(direction, n1, n2, n3, scale, 
		  (fdr_complex *)data, ldx, ldx2, 
		  (fdr_complex *)data, ldy, ldy2,
          table, work, isys);

  if(verbose>1) fprintf(stdout,"Completed 3DFFT\n");

  if(use_fftshift)
	cfftshift(n3,n2,n1,data);

  if(verbose>1) fprintf(stdout, "Completed FFTSHIFT\n");

  /* Free the temporary memory */
  free(table); free(work);

  if(verbose>1) fprintf(stdout, "**************Returning from 3DFFT*****************\n");
  return 1;

}

int psfstack_fft(long n3, long n2, long n1, fdr_complex *data) 
{
  int result;
  long ldx, ldy, i;
  long maxdim;

  fdr_complex *tmpdata = NULL;
  float *table = NULL;
  float *work = NULL;
  int isys[2];

  if(data == NULL) {
	fprintf(stderr, "psfstack_fft called with NULL data pointers.\n");
	return 0;
  }
  
  /* Set up the dimensions */
  maxdim = n1>n2 ? n1:n2;
  ldx = n1;
  ldy = n1;

  /* Allocate memory */
  table = (float *)calloc(((2*n1+256) + (2*n2+256)),sizeof(float));
  work = (float *)calloc((2*maxdim),sizeof(float));
  if(table == NULL || work == NULL) {
	fprintf(stderr, "Could not allocate memory!\n");
	return 0;
  }
  if(verbose>1) fprintf(stdout,"Allocated memory\n");

  tmpdata = data;
  for(i=0;i<n3;i++,tmpdata+=n2*ldx) {

	//fprintf(stdout, "As a test ... this first voxel is %f+%fi\n", tmpdata[0].re, tmpdata[0].im);
	if(verbose>1) fprintf(stdout,"Performing FFT of slice %d of %d\n",i,n3);
	if(use_fftshift)
	  c2dfftshift(n2,n1,tmpdata);
	/* Do the 3D IFFT using SCSL */
	isys[0] = 1;
	ccfft2d(0, n1, n2, 1.0f, 
			(fdr_complex *)tmpdata, ldx, 
			(fdr_complex *)tmpdata, ldy,
			table, work, isys);
	// if(verbose>1) fprintf(stdout, "Prepped table\n");
	ccfft2d(1, n1, n2, 1.0f, 
			(fdr_complex *)tmpdata, ldx, 
			(fdr_complex *)tmpdata, ldy,
			table, work, isys);
	if(use_fftshift)
	  c2dfftshift(n2,n1,tmpdata);

  }
  if(verbose>1) fprintf(stdout,"Completed 3DFFT\n");

  /* Free the memory we've been using */
  free(table);free(work);

  /* Finally, we're successful */
  if(verbose>1) fprintf(stdout, "************ Returning from psfstack_fft ***********\n");
  return 1;

}

int normalize_complexpsfstack(long n3, long n2, long n1, fdr_complex *data) 
{
  
  /* NOTE: 
   * This program assumes that you've got no information in imaginary space
   * in other words, this is just a memory saving feature!
   * real-valued data only!
   */
  int result;
  long i, j, k;
  long index;
  float tempsum, maxsum;

  if(data == NULL) {
	fprintf(stderr, "psfstack_fft called with NULL data pointers.\n");
	return 0;
  }

  maxsum = 0;
  /* Iterate over depth point */
  for(k=0;k<n3;k++) {
	tempsum = 0;
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		tempsum += data[index].re;
	  }
	}
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		data[index].re /= tempsum;
	  }
	}
	if(tempsum>maxsum) maxsum = tempsum;
  }

  fprintf(stdout,"NORMALIZE: Maximum sum is: %f\n", maxsum);

  /* I can't figure out why this was ever the case - Mar 17 */
  /* Now we know the max, let's normalize */
  //#pragma omp parallel for
  //for(i=0;i<n3*n2*n1;i++) {
  //	data[i] /= maxsum;
  //}


  if(verbose>1) fprintf(stdout,"Completed PSF stack normalization\n");

  /* Finally, we're successful */
  return 1;

}


int normalize_psfstack(long n3, long n2, long n1, float *data) 
{
  // Changed Mar 17 to normalize each slice to 1
  int result;
  long i, j, k;
  long index;
  float tempsum, maxsum;

  if(data == NULL) {
	fprintf(stderr, "psfstack_fft called with NULL data pointers.\n");
	return 0;
  }

  maxsum = 0;
  /* Iterate over depth point */
  for(k=0;k<n3;k++) {
	tempsum = 0;
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		tempsum += data[index];
	  }
	}
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		data[index] /= tempsum;
	  }
	}
	if(tempsum>maxsum) maxsum = tempsum;
  }

  fprintf(stdout,"NORMALIZE: Maximum sum is: %f\n", maxsum);

  /* I can't figure out why this was ever the case - Mar 17 */
  /* Now we know the max, let's normalize */
  //#pragma omp parallel for
  //for(i=0;i<n3*n2*n1;i++) {
  //	data[i] /= maxsum;
  //}


  if(verbose>1) fprintf(stdout,"Completed PSF stack normalization\n");

  /* Finally, we're successful */
  return 1;

}

int limit_recovery_filter(long nelems, float limit, fdr_complex *filterdata) {

  long i;
  float val1, val2; //, angle;
  

  if(filterdata == NULL) {
	fprintf(stderr, "limit_recovery_filter called with NULL data\n");
	return 0;
  }

#pragma omp parallel for private(val1,val2)
  for(i=0;i<nelems;i++) {
	val1 = sqrt(filterdata[i].re*filterdata[i].re + filterdata[i].im*filterdata[i].im);
	// angle = arctan(filterdata[i].im/filterdata[i].re);
	/*
	if(val1 > limit) {
	  val2 = limit/val1;
	  filterdata[i].re = filterdata[i].re * val2;
	  filterdata[i].im = filterdata[i].im * val2;
	}
	*/
	
	if(val1 > limit/2) {
	  val2 = limit - limit/2*exp(-(val1-limit/2)/limit);
	  filterdata[i].re = filterdata[i].re * val2/val1;
	  filterdata[i].im = filterdata[i].im * val2/val1;
	}
	
  }
  
  if(verbose>1) fprintf(stdout, "******** Completed limit_recovery_filter *********\n");
  return 1;

}


/***************************************************************
 ***************************************************************
 ******************* DATA COMPARE FUNCTIONS ********************
 ***************************************************************
 ***************************************************************/


int compare_complex_data(long n3,long n2,long n1, fdr_complex *data1, fdr_complex *data2) {

  long j,k,i, index;
  
  fprintf(stdout, "Beginning comparison...\n");
  for(k=0;k<n3;k++) {
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		if(fabs(data1[index].re) > fabs(data2[index].re))
		  fprintf(stdout, "At k=%d,j=%d,i=%d, REAL data1 is larger than data2 (%f vs %f)\n", k,j,i,data1[index].re, data2[index].re);
		if(fabs(data1[index].im) > fabs(data2[index].im))
		  fprintf(stdout, "At k=%d,j=%d,i=%d, REAL data1 is larger than data2 (%f vs %f)\n", k,j,i,data1[index].im, data2[index].im);
	  }
	}
  }
  fprintf(stdout, "Ending comparison.  If you don't see any errors, you're good!\n");

  return 1;

}

int compare_real_data(long n3,long n2,long n1, float *data1, float *data2) {

  long j,k,i, index;
  
  fprintf(stdout, "Beginning comparison...\n");
  for(k=0;k<n3;k++) {
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		if(data1[index] > data2[index]) {
		  fprintf(stdout, "At k=%d,j=%d,i=%d, data1 is larger than data2 (%f vs %f)\n", k,j,i,data1[index], data2[index]);
		}
	  }
	}
  }
  fprintf(stdout, "Ending comparison.  If you don't see any errors, you're good!\n");

  return 1;

}

/***************************************************************
 ***************************************************************
 ******************* FILTER CONSTRUCTION ***********************
 ***************************************************************
 ***************************************************************/


int build_rolloff_filter(long n3, long n2, long n1, float weight, float maxslope, float *outdata) {

  long i,j,k,index1,index2;
  float startk,starti,curk,curi,dk,di,slope;
  float *rolloff;
  
  if(outdata == NULL) {
	fprintf(stderr, "build_rolloff_filter called with NULL data\n");
	return 0;
  }

  rolloff = calloc(n3*n1,sizeof(float));
  if(rolloff == NULL) {
	fprintf(stderr, "Could not allocate memory to build rolloff filter\n");
	return 0;
  }

  /* Do some calculations outside the loop */
  // Changed July 6 to account for even and odd length data sets
  dk = 1/(2*M_PI);
  di = 1/(2.0);
  /*
  startk = (-(n3-1)/2.0) / (2.0*M_PI) - dk/2.;
  starti = (-(n1-1)/2.0) / (2.0) - di/2.;
  */
  if(n3%2)
    startk = (-(n3-1.0)/2.0) / (2.0*M_PI);
  else
    startk = (-(  n3  )/2.0) / (2.0*M_PI);
  if(n1%2)
 		starti = (-(n1-1.0)/2.0) / (2.0);
 	else
 		starti = (-(  n1  )/2.0) / (2.0);
	  

  for(k=0,curk=startk;k<n3;k++,curk+=dk) {
	for(i=0,curi=starti;i<n1;i++,curi+=di) {
	  index1 = k*n1 + i;
	  slope = -curk/curi;
	  /* Special case -- along the zero-axis */
	  //if( (fabs(curk)<=1e-8) || fabs(curi)<=(1e-8)) {
	  if(fabs(curk)<=(1e-8)) {
		rolloff[index1] = 1.0;
		continue;
	  }
	  if(fabs(curi)<=(1e-8) && fabs(curk) < 0.5) {
		rolloff[index1] = 1.0;
		continue;
	  }
	  /* Should be no data in these areas */
	  if(fabs(slope) >= 1.0) {
		rolloff[index1] = 0.0;
		continue;
	  }
	  /* Cosine rolloff from maxslope down to 1.0 */
	  if(slope >= maxslope) {
		rolloff[index1] =pow(cos(M_PI/2 * (fabs(slope)-maxslope) / (1.0-maxslope)), 2);
		continue;
	  }
	  /* Zero out the out of focus data beyond weight */
	  if(slope <= -weight) {
		rolloff[index1] = 0.0;
		continue;
	  }
	  /* Want to keep all data in the good quadrant*/
	  if(slope >= 0) {
		rolloff[index1] = 1.0;
		continue;
	  }
	  /* Everything else should be eliminated. We should be in rolloff region. */
	  /* It is 1.0 at slope = 0.0, and 0.0 at slope = weight */
	  rolloff[index1] = pow(cos( M_PI/2.0 * slope/weight ),2);
	  /*
	  if(slope>0)
		rolloff[index1] = 1.0;
	  else if (slope<-1.0)
		rolloff[index1] = 1.0;  // change this to zero if you want that whole quadrant gone
	  else if(slope<=-weight)
		rolloff[index1] = 0.0; 
	  else {
		rolloff[index1] = pow(cos( M_PI/2.0 * slope/weight ),2);
		// fprintf(stdout, "We are setting %d,%d to be %f (from %f)\n", k,i,rolloff[index1],slope);
	  }
	  */
	}

  }

  for(k=0;k<n3;k++) {
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index1 = k*n2*n1 + j*n1 + i;
		index2 = k*n1 + i;
		outdata[index1] = rolloff[index2];
	  }
	}
  }

  free(rolloff);

  if(verbose>1) fprintf(stdout, "****** Returning from build_rolloff_filter ******\n");

  return 1;

}

int build_wiener_filter(long n3, long n2, long n1, long newelems, fdr_complex *fftdata, float *outdata)
{
  long nelems,i,j,k, index;
  int result;
  float tmp, noise,noise2;
  float *ptmp;
  float *rad, *count;
  float radius, frac;
  float z,y,x;
  float scale;
  long iradius;

  if(fftdata == NULL) {
	fprintf(stderr, "build_noise_dampening_filter called with NULL data\n");
	return 0;
  }

  rad = calloc(newelems,sizeof(float));
  count = calloc(newelems,sizeof(float));
  if(rad == NULL || count == NULL) {
	fprintf(stderr, "Could not allocate memory for radial power spectrum\n");
	return 0;
  }

  nelems = n3*n2*n1;

  /* Calculate the power spectrum */
#pragma omp parallel for
  for(i=0;i<nelems;i++) {
	outdata[i] = (float)(fftdata[i].re*fftdata[i].re + fftdata[i].im*fftdata[i].im);
  }
  if(verbose>1) fprintf(stdout, "Calculated power spectrum.\n");

  /* Let's get a mean value for the noise floor */
  /* For this we'll use a bit of data just off the center line */
  noise = 0;
  for(k=n3/2-50;k<n3/2-47;k++) {
	for(j=0;j<3;j++) {
	  for(i=0;i<3;i++) {
		index = k*n2*n1 + j*n2 + i;
		// fprintf(stdout, "Index is: %d\n", index);
		noise += outdata[index];
	  }
	}
  }
  noise /= (3.*3.*3.);
  if(verbose>1) fprintf(stdout, "Calculated noise floor of %f\n", noise);

  scale = (newelems-2) / sqrt( (n3-1.)/2.*(n3-1)/2. + (n2-1.)/2.*(n2-1)/2. + (n1-1.)/2.*(n1-1.)/2. );
  /* Now average over the radius to average out noise */
#pragma omp parallel for private (z,y,x,index,radius,iradius,frac,j,i)
  for(k=0;k<n3;k++) {
  	if(n3%2)
			z = -(n3-1)/2. + k;
		else
			z = -( n3 )/2. + k;
		for(j=0;j<n2;j++) {
			if(n2%2)
			  y = -(n2-1)/2. + j;
			else
				y = -( n2 )/2. + j;
		  for(i=0;i<n1;i++) {
		  	if(n1%2)
					x = -(n1-1)/2. + i;
				else
				  x = -( n1 )/2. + i;
				index = k*n2*n1 + j*n1 + i;
				radius = sqrt(z*z + y*y + x*x) * scale;
				iradius = (int)floor(radius);
				frac = radius - iradius;
				if(iradius < newelems-1) {
				  rad[iradius] += (1-frac)*outdata[index];
				  count[iradius] += (1-frac);
				}
				if(iradius+1 < newelems-1) {
				  rad[iradius+1] += frac*outdata[index];
				  count[iradius+1] += frac;
				}
		  }
		}
  }

  noise2 = 0;
  for(i=0;i<9;i++) {
	noise2 += rad[newelems-i];
  }
  noise2 /= 9;

  if(verbose>1) fprintf(stdout, "Calculated alternate noise floor of %f\n", noise2);

#pragma omp parallel for
  for(i=0;i<newelems;i++) {

	//fprintf(stdout, "PS At %d Rad is %f and count is %f, resulting in %f\n", i,rad[i], count[i],rad[i]/count[i]);

	if(count[i] == 0) 
	  rad[i] = 0;
	else
	  rad[i] = rad[i]/count[i];
  }

  /* Now we can build the filter itself */
  /* This is the radial filter, have to copy it to 3d later */
#pragma omp parallel for
  for(i=0;i<newelems;i++) {
	/* This is the traditional Wiener filter */
	/* Our "signal" is just power spec of image minus noise */
	if(rad[i] != 0) {
	  rad[i] = 1. /  (1. + noise / (rad[i]-noise) );
	  if(rad[i] < 0)
		rad[i] = 0;
	}
	//fprintf(stdout, "WF At %d filter is %f\n", i, rad[i]);
  }

#pragma omp parallel for private (z,y,x,index,radius,iradius,frac,j,i)
  for(k=0;k<n3;k++) {
  	if(n3%2)
			z = -(n3-1)/2. + k;
		else
			z = -( n3 )/2. + k;
		for(j=0;j<n2;j++) {
			if(n2%2)
			  y = -(n2-1)/2. + j;
			else
				y = -( n2 )/2. + j;
		  for(i=0;i<n1;i++) {
		  	if(n1%2)
					x = -(n1-1)/2. + i;
				else
					x = -( n1 )/2. + i;
				index = k*n2*n1 + j*n1 + i;
				radius = sqrt(z*z + y*y + x*x) * scale;
				iradius = (int)floor(radius);
				frac = radius - iradius;
				outdata[index] = (1-frac)*rad[iradius] + frac*rad[iradius+1];
				if(outdata[index] < 0) {
				  fprintf(stdout,"Bug?  frac=%f r=%f ir=%d rad=%f rad2=%f\n", frac, radius, iradius, rad[iradius],rad[iradius+1]);
				  outdata[index] = 0;
				}
				if(outdata[index] > 1.0) {
				  fprintf(stdout, "Bug?  frac=%f r=%f ir=%d rad=%f rad2=%f\n", frac, radius, iradius, rad[iradius],rad[iradius+1]);
				  outdata[index] = 1;
				}
		  }
		}
  }

  if(verbose>1) fprintf(stdout, "Calculated wiener filter\n");

  return 1;

}

int build_antialiased_data(long n3, long n2, long n1, long cutoff, fdr_complex *fftdata, fdr_complex *outdata) 
{
  long i, j, k, nv, t, index1, index2;
  int result;
  
  if(fftdata == NULL || outdata == NULL) {
	fprintf(stderr, "build_antialiased_data called with NULL data pointers.\n");
	return 0;
  }


  // views*3
  for(t=0;t<3;t++) {
		for(k=0;k<n3;k++) {
		  nv = k + n3*t;
		  if(nv < 3*n3/2 - cutoff/2) continue;
		  if(nv > 3*n3/2 + cutoff/2) continue;
		  
		  fprintf(stdout, "Remapping k=%d to nv=%d\n", k, nv);
		  // slices
		  for(j=0;j<n2;j++) {
				// detectors
				for(i=0;i<n1;i++) {
				  index1 = nv*n2*n1 + j*n1 + i;
				  index2 = k*n2*n1 + j*n1 + i;
				  outdata[index1].re = fftdata[index2].re;
				  outdata[index1].im = fftdata[index2].im;
				  //(fdr_complex)outdata[index1] = (fdr_complex)fftdata[index2];
				  /*
					if(random()) {
					fprintf(stdout, "Setting %d,%d,%d from %d,%d,%d\n", k,j,i, nv,j,i);
					fprintf(stdout, "This value is %f,%f from %f,%f\n", outdata[index1].re, 
					outdata[index1].im, fftdata[index2].re,fftdata[index2].im);
					}
				  */
				}
		  }
		}
  }

  return 1;
  
}

int build_noise_dampening_filter(long n3, long n2, long n1, fdr_complex *fftdata, float *outdata)
{

  long nelems,i,j,k, index;
  int result;
  float tmp, noise;
  float *ptmp;
  float mx;

  mx = 0;

  if(fftdata == NULL) {
	fprintf(stderr, "build_noise_dampening_filter called with NULL data\n");
	return 0;
  }
  
  nelems = n3*n2*n1;

  /* Calculate the power spectrum */
#pragma omp parallel for
  for(i=0;i<nelems;i++) {
	outdata[i] = (float)(fftdata[i].re*fftdata[i].re + fftdata[i].im*fftdata[i].im);
  }
  if(verbose>1) fprintf(stdout, "Calculated power spectrum.\n");

  /* Let's get a mean value for the noise floor */
  noise = 0;
  for(k=0;k<10;k++) {
	for(j=0;j<10;j++) {
	  for(i=0;i<10;i++) {
		index = k*10*10 + j*10 + i;
		noise += outdata[index];
	  }
	}
  }
  noise /= (10.*10.*10.);
  if(verbose>1) fprintf(stdout, "Calculated noise floor of %f\n", noise);

  /* Now we can build the filter itself */
#pragma omp parallel for private(mx)
  for(i=0;i<nelems;i++) {
	/* This is the traditional Wiener filter */
	outdata[i] = (outdata[i]-noise) /  ( (outdata[i]-noise)+noise/(outdata[i]-noise) );
	/* This is the RMH suggested noise dampening filter */
	// outdata[i] = outdata[i]/noise-1.0;
	mx = outdata[i] > mx ? outdata[i] : mx;
  }

#pragma omp parallel for
  for(i=0;i<nelems;i++)
	outdata[i] = outdata[i] / mx;

  if(verbose>1) fprintf(stdout, "Calculated noise dampening filter\n");

  return 1;

}

int build_bw_filter(long n3, long n2, long n1, float bandlimit, float *bwfilterdata)
{

  long k,j,i, index;
  float startk,dk, startj,dj, starti,di, curk, curj, curi;
  float bw_i, bw_j, scale;

  if(bwfilterdata == NULL) {
	fprintf(stderr, "build_bw_filter called with NULL data pointers.\n");
	return 0;
  }

	// Changed July 6 to take odd and even numbered bins into account
  //startk = (-(n3-1.0)/2.0) / (2.0*M_PI);
  if(n3%2)
  	startk = (-(n3-1.0)/2.0) / (2.0*M_PI);
  else
  	startk = (-(  n3  )/2.0) / (2.0*M_PI);
  dk = 1/(2*M_PI);
  //startj = (-(n2-1.0)/2.0) / (2.0);
  if(n2%2)
  	startj = (-(n2-1.0)/2.0) / (2.0);
  else
  	startj = (-(  n2  )/2.0) / (2.0);
  dj = 1/(2.0);
  //starti = (-(n1-1.0)/2.0) / (2.0);
  if(n1%2)
  	starti = (-(n1-1.0)/2.0) / (2.0);
  else
  	starti = (-(  n1  )/2.0) / (2.0);
  di = 1/(2.0);

  for(k=0,curk=startk;k<n3;k++,curk+=dk) {
		for(j=0,curj=startj;j<n2;j++,curj+=dj) {
	  	for(i=0,curi=starti;i<n1;i++,curi+=di) {
			index = k*n2*n1 + j*n1 + i;
		
			// At this point we want to check the bandlimit of the data
			// Three cases: 
			//     less than bandlimit*0.9 ==> scale=1
			//     less than bandlimit     ==> scale=cosine rolloff
			//     greater than bandlimit  ==> scale = 0
			bw_i = bandlimit*(n1-1)/4.0;
			bw_j = bandlimit*(n2-1)/4.0;
			// if(verbose>1) fprintf(stdout, "Bandlimit i is %f and j is %f at %f and %f\n", bw_i, bw_j, curi, curj);
			
			if(fabs(curi) > bw_i || fabs(curj) > bw_j) {
			  scale = 0.0;
			}
			else if(fabs(curi) < 0.9*bw_i && fabs(curj) < 0.9*bw_j) {
			  scale = 1.0;
			}
			
			else if( (fabs(curi) >= 0.9*bw_i && fabs(curi) <= bw_i) && 
					 (fabs(curj) >= 0.9*bw_j && fabs(curj) <= bw_j) ) {
			  scale = pow(cos(M_PI/2 * (fabs(curi)-0.9*bw_i)/(0.1*bw_i)),2) * 
				pow(cos(M_PI/2 * (fabs(curj)-0.9*bw_j)/(0.1*bw_j)),2);
			}
			
			else if(fabs(curi) >= 0.9*bw_i && fabs(curi) <= bw_i) {
			  //tmp = 1-(bw_i- curi)/(0.1* bw_i);
			  //scale = pow(cos(tmp*2*M_PI/4.0),2);
			  scale = pow(cos(M_PI/2 * (fabs(curi)-0.9*bw_i)/(0.1*bw_i)),2);
			}
			else if(fabs(curj) >= 0.9*bw_j && fabs(curj) <= bw_j) {
			  //tmp = 1-(bw_j -curj)/(0.1*bw_j);
			  //scale = pow(cos(tmp*2*M_PI/4.0),2);
			  scale = pow(cos(M_PI/2 * (fabs(curj)-0.9*bw_j)/(0.1*bw_j)),2);
			}
			bwfilterdata[index] = scale;
		  }
		}
  }

  if(verbose) fprintf(stdout, "***** Completed build_bw_filter *****\n");

}

int resample_for_comparison(long psfn3, long psfn2, long psfn1, float *psfcomparedata,
                            long    n3, long    n2, long    n1, fdr_complex *sgfftdata)
{
  /* 
   * In this function we have the opposite case of build_fdp - we have the sinogram fft data,
   * and we wish to resample it to look like the stack of 2D PSF FFTs. 
   *
   */

  long nelems;
  long psfindex;
  long k,j,i, index1,index2;
  float startk,dk, startj,dj, starti,di, curk, curj, curi;
  float mindepth, maxdepth, ddepth, frac;
  long idepth;
  float depth, slope;
  float tmp1, tmp2, denom1, denom2;
  float bw_i, bw_j, scale, tmp, abs1, abs2;
  float minslope, maxslope, dslope, minview, maxview, dview, mindet, maxdet, ddet;
  float curdet, curslope, z_pos, iz;

  if(psfcomparedata == NULL || sgfftdata == NULL) {
	fprintf(stderr, "resample_for_comparison called with NULL data pointers.\n");
	return 0;
  }

  nelems = n3*n2*n1;

  /* For the stack of PSF FFTs:
   * It is assumed that z=0 is the slope=-1 depth point 
   * and that z=nx is the slope=+1 depth point 
   */
 
  /* Let's make sure that the data is appropriate before we go any further */
  /* The X and Y dimensions should match up */
  if(psfn2 != n2 || psfn1 != n1) {
	fprintf(stderr, "The X,Y dimensions of the point spread function and the sinogram data do not match.\n");
	fprintf(stderr, "PSF: Z=%d Y=%d X=%d, SG: Z=%d, Y=%d, X=%d\n", psfn3, psfn2, psfn1, n3, n2, n1);
	return 0;
  }

  /* So now we come to the guts of the operation
   * At this point we are taking the 3d analog to the bow-tie and turning it into
   * data to compare to the 2D stack of PSF FFTs
   *
   * Coordinates are Z=View, Y=Slice, X=Detector
   * In Fourier Space, we are looking at creating a filter H(Rx,Ry,-PHI/Rx)
   * so Coordinates are:
   *
   * Z=-PHI/Rx, Y=Ry, X=-PHI/Rx
   *
   * This puts the calculation of the slope in the innermost bracket
   * This could like be optimized by doing those calculations first?
   */

  /* Do some calculations outside the loop */
  minslope = -1.0;
  maxslope = 1.0;
  dslope = (maxslope - minslope)/(psfn3-1);

  /*
  minview = (-(n3-1.0)) / (2.0*M_PI); 
  maxview = ((n3-1.0)) / (2.0*M_PI);
  dview = (maxview-minview)/(n3);

  mindet = (-(n1-1.0)) / (2.0);
  maxdet = ((n1-1.0)) / (2.0);
  ddet = (maxdet - mindet)/(n1);
  */
  if(n2%2)
    minview = -(n3-1)/2.0 / (2*M_PI);
  else
    minview = - n3 / 2.0 / (2*M_PI);
  dview = 1/(2*M_PI);
  maxview = minview + n3*dview;
  
  if(n1%2)
    mindet = -(n1-1)/2.0/2.0;
  else
    mindet = -n1 / 2.0 / 2.0;
  ddet = 1/2.0;
  maxdet = mindet + n1*ddet;
  
  /*
  startk = (-(n3-1.0)/2.0) / (2.0*M_PI);
  dk = 1/(2*M_PI);
  startj = (-(n2-1.0)/2.0) / (2.0);
  dj = 1/(2.0);
  starti = (-(n1-1.0)/2.0) / (2.0);
  di = 1/(2.0);
  ddepth = 2.0/ (psfn3-1);
  */
  maxdepth = 1.0;
  mindepth = -1.0;
  fprintf(stdout, "Details are: startk %f dk %f startj %f dk %f starti %f di %f\n", startk, dk, startj,dj,starti,di);
  fprintf(stdout, "Depth details are: max %f min %f ddepth %f\n", maxdepth, mindepth, ddepth);
  fprintf(stdout, "slope: %f : %f : %f  views: %f : %f : %f  dets: %f : %f : %f \n", 
          minslope, dslope, maxslope, minview, dview, maxview, mindet, ddet, maxdet);
  //fprintf(stdout, "Bandlimit at: %f\n", bandlimit);
  fprintf(stdout, "PSF DIMS: Z=%d Y=%d X=%d\n", psfn3, psfn2, psfn1);
  fprintf(stdout, "SG  DIMS: Z=%d Y=%d X=%d\n", n3, n2, n1);
  /*
  for(i=0,curdet=mindet;i<n1;i++,curdet+=ddet) {
    fprintf(stdout, "at i=%d, det=%f\n", i, curdet);
  }
	*/
  /* Loop over all elements to get the appropriate information */
  /* Note this is iterating in frequency space, so k,j,i represent the coordinates denoted above */
  for(k=0,curslope=minslope;k<psfn3;k++,curslope+=dslope) {
	if(verbose) fprintf(stdout, "Now sampling Z=%d, slope=%f\n", k, curslope);
	for(j=0;j<n2;j++) {
 	  for(i=0,curdet=mindet;i<n1;i++,curdet+=ddet) {
        psfindex =  k*psfn2*psfn1 + j*psfn1 + i;
        z_pos = -curslope * curdet;
        iz = (z_pos - minview)/dview;
        frac = iz - floor(iz);
				index1 =   (int)iz  *n2*n1 + j*n1 + i;
        index2 = ((int)iz+1)*n2*n1 + j*n1 + i;
        
		if(index1>nelems || index2 > nelems || index1 < 0 || index2 < 0) {
		  //fprintf(stdout, "WHOA!  index > nelems, this shouldn't be possible! %d %d, from %d %d %d\n", index1, index2, (int)iz, j, i);
		  continue;
		}

		/*
		  fprintf(stdout, "At coords Z=%d Y=%d X=%d we have curk=%f curj=%f, curi=%f\n", k, j, i, curk, curj, curi);
        */
        //fprintf(stdout, "At coords Z=%d Y=%d X=%d we have curslope=%f curdet=%f, zpos=%f\n", k, j, i, curslope, curdet, z_pos);
          /*
		  fprintf(stdout, "Slope is %f at Z=%d, Y=%d, X=%d\n", slope, k, j, i);
		  fprintf(stdout, "This corresponds to a depth of %f (index %d)\n", depth, idepth);
		  fprintf(stdout, "And a psfindex of %d and %d (max of %d)\n", psfindex1, psfindex2, psfn3*psfn2*psfn1);
        */
		/* First, let's interpolate to get the right frequency space value */
        
        abs1 = sqrt( sgfftdata[index1].re*sgfftdata[index1].re + sgfftdata[index1].im*sgfftdata[index1].im);
        abs2 = sqrt( sgfftdata[index2].re*sgfftdata[index2].re + sgfftdata[index2].im*sgfftdata[index2].im);
        psfcomparedata[psfindex] += (1-frac)*abs1 + frac*abs2;
		//tmp1 = (1-frac)*psffftdata[psfindex1].re + frac*psffftdata[psfindex2].re;
		//tmp2 = (1-frac)*psffftdata[psfindex1].im + frac*psffftdata[psfindex2].im;
		
        //fdpfilterdata[index].re = tmp1 * scale;
        //fdpfilterdata[index].im = tmp2 * scale;
	  }
	}
  }

  if(verbose) fprintf(stdout,"Completed building the compare to PSF filter.\n");
  
  /* Finally, we're successful */
  if(verbose) fprintf(stdout, "*********** Returning from reample_for_comparison **************\n");
  
  return 1;


}

int build_fdp(long psfn3, long psfn2, long psfn1, fdr_complex *psffftdata, 
			  long    n3, long    n2, long    n1, fdr_complex *fdpfilterdata, 
			  float mxsl, float bandlimit, int flags)
{

  long nelems;
  long psfindex1,psfindex2;
  long k,j,i, index;
  float startk,dk, startj,dj, starti,di, curk, curj, curi;
  float mindepth, maxdepth, ddepth, frac;
  long idepth;
  float depth, slope;
  float tmp1, tmp2, denom1, denom2;
  float bw_i, bw_j, scale, tmp;

  fprintf(stdout, "Maxslope is %f\n", mxsl);

  if(psffftdata == NULL || fdpfilterdata == NULL) {
	fprintf(stderr, "build_fdp called with NULL data pointers.\n");
	return 0;
  }

  nelems = n3*n2*n1;

  /* For the stack of PSF FFTs:
   * It is assumed that z=0 is the slope=-1 depth point 
   * and that z=nx is the slope=+1 depth point 
   */
 
  /* Let's make sure that the data is appropriate before we go any further */
  /* The X and Y dimensions should match up */
  if(psfn2 != n2 || psfn1 != n1) {
	fprintf(stderr, "The X,Y dimensions of the point spread function and the sinogram data do not match.\n");
	fprintf(stderr, "PSF: Z=%d Y=%d X=%d, SG: Z=%d, Y=%d, X=%d\n", psfn3, psfn2, psfn1, n3, n2, n1);
	return 0;
  }

  /* So now we come to the guts of the operation
   * At this point we are taking the stack of PSF FFT's and transforming
   * it into the depth-dependent point spread function - the 3d analog to the bowtie - 
   * that will then be used to deconvolve the depth dependent point spread function
   * of the lens
   *
   * Coordinates are Z=View, Y=Slice, X=Detector
   * In Fourier Space, we are looking at creating a filter H(Rx,Ry,-PHI/Rx)
   * so Coordinates are:
   *
   * Z=-PHI/Rx, Y=Ry, X=-PHI/Rx
   *
   * This puts the calculation of the slope in the innermost bracket
   * This could like be optimized by doing those calculations first?
   */

  /* Do some calculations outside the loop */
  // Changed July 6 2006 to account for odd and even numbered bins
  //startk = (-(n3-1.0)/2.0) / (2.0*M_PI);
  if(n3%2)
  	startk = (-(n3-1.0)/2.0) / (2.0*M_PI);
  else
  	startk = (-(  n3  )/2.0) / (2.0*M_PI);
  dk = 1/(2*M_PI);
  //startj = (-(n2-1.0)/2.0) / (2.0);
  if(n2%2)
  	startj = (-(n2-1.0)/2.0) / (2.0);
  else
  	startj = (-(  n2  )/2.0) / (2.0);
  dj = 1/(2.0);
  //starti = (-(n1-1.0)/2.0) / (2.0);
  if(n1%2)
  	starti = (-(n1-1.0)/2.0) / (2.0);
  else
  	starti = (-(  n1  )/2.0) / (2.0);
  di = 1/(2.0);
  
  ddepth = 2.0/ (psfn3-1);
  maxdepth = mxsl;
  mindepth = -mxsl;
  fprintf(stdout, "Details are: startk %f dk %f startj %f dk %f starti %f di %f\n", startk, dk, startj,dj,starti,di);
  fprintf(stdout, "Depth details are: max %f min %f ddepth %f\n", maxdepth, mindepth, ddepth);
  fprintf(stdout, "Bandlimit at: %f\n", bandlimit);
  fprintf(stdout, "PSF DIMS: Z=%d Y=%d X=%d\n", psfn3, psfn2, psfn1);
  fprintf(stdout, "SG  DIMS: Z=%d Y=%d X=%d\n", n3, n2, n1);
  

  /* Loop over all elements to get the appropriate information */
  /* Note this is iterating in frequency space, so k,j,i represent the coordinates denoted above */
  for(k=0,curk=startk;k<n3;k++,curk+=dk) {
	if(verbose) fprintf(stdout, "Now sampling Z=%d\n", k);
	for(j=0,curj=startj;j<n2;j++,curj+=dj) {
	  for(i=0,curi=starti;i<n1;i++,curi+=di) {
		index = k*n2*n1 + j*n1 + i;
		if(index>nelems) {
		  fprintf(stdout, "WHOA!  index > nelems, this shouldn't be possible!\n");
		  continue;
		}
		slope = -curk/curi;
		scale = 1.0;
		/* Just zero out everything that isn't "interesting */
		
		if(fabs(slope) >= 1.0) {
		  fdpfilterdata[index].re = 1.0;
		  fdpfilterdata[index].im = 0.0;
		  continue;
		}
		
		
		/*
		if(slope > maxdepth || slope < mindepth) {
		  if(slope < 1.0 || slope > -1.0) 
		    fdpfilterdata[index].re = 1.0; // change this to zero if you want that whole anti-bowtie gone
		  else
			fdpfilterdata[index].re = 0.0;
		  fdpfilterdata[index].im = 0.0;
		  continue;
		}
		*/
		
		

		// Cosine roll off from maxslope down to 1.0
		/*
		if(fabs(slope) >= maxdepth) {
		  scale = pow(cos(M_PI/2 * (fabs(slope)-maxdepth) / (1.0-maxdepth)), 2);
		}
		*/
		

		
		/*
		// At this point we want to check the bandlimit of the data
		// Three cases: 
		//     less than bandlimit*0.9 ==> scale=1
		//     less than bandlimit     ==> scale=cosine rolloff
		//     greater than bandlimit  ==> scale = 0
		bw_i = bandlimit*(n1-1)/4.0;
		bw_j = bandlimit*(n2-1)/4.0;
		//fprintf(stdout, "Bandlimit i is %f and j is %f at %f and %f\n", bw_i, bw_j, curi, curj);
		
		if(fabs(curi) > bw_i || fabs(curj) > bw_j) {
		  scale = 0.0;
		}
		else if(fabs(curi) < 0.9*bw_i && fabs(curj) < 0.9*bw_j) {
		  scale = 1.0;
		}
		
		else if( (fabs(curi) >= 0.9*bw_i && fabs(curi) <= bw_i) && 
				 (fabs(curj) >= 0.9*bw_j && fabs(curj) <= bw_j) ) {
		  scale = pow(cos(M_PI/2 * (fabs(curi)-0.9*bw_i)/(0.1*bw_i)),2) * 
			pow(cos(M_PI/2 * (fabs(curj)-0.9*bw_j)/(0.1*bw_j)),2);
		}
		
		else if(fabs(curi) >= 0.9*bw_i && fabs(curi) <= bw_i) {
		  //tmp = 1-(bw_i- curi)/(0.1* bw_i);
		  //scale = pow(cos(tmp*2*M_PI/4.0),2);
	 	  scale = pow(cos(M_PI/2 * (fabs(curi)-0.9*bw_i)/(0.1*bw_i)),2);
 		}
		else if(fabs(curj) >= 0.9*bw_j && fabs(curj) <= bw_j) {
		  //tmp = 1-(bw_j -curj)/(0.1*bw_j);
		  //scale = pow(cos(tmp*2*M_PI/4.0),2);
		  scale = pow(cos(M_PI/2 * (fabs(curj)-0.9*bw_j)/(0.1*bw_j)),2);
	   }
		*/
		// depth = (slope-mindepth);  // what the fuck is this line in here for?  idiot
		depth = slope + 1.0;
		idepth = floor(depth/ddepth); 
		frac = depth/ddepth-idepth;
		psfindex1 =  idepth   *psfn2*psfn1 + j*psfn1 + i;
		psfindex2 = (idepth+1)*psfn2*psfn1 + j*psfn1 + i;
		/*
		  fprintf(stdout, "At coords Z=%d Y=%d X=%d we have curk=%f curj=%f, curi=%f\n", k, j, i, curk, curj, curi);
		  
		  fprintf(stdout, "Slope is %f at Z=%d, Y=%d, X=%d\n", slope, k, j, i);
		  fprintf(stdout, "This corresponds to a depth of %f (index %d)\n", depth, idepth);
		  fprintf(stdout, "And a psfindex of %d and %d\n", psfindex1, psfindex2);
		*/
		/* First, let's interpolate to get the right frequency space value */
		tmp1 = (1-frac)*psffftdata[psfindex1].re + frac*psffftdata[psfindex2].re;
		tmp2 = (1-frac)*psffftdata[psfindex1].im + frac*psffftdata[psfindex2].im;
		
		if(flags == INVERSE_FILTER) {
		  /* 1 / (a+bi) = a/(a^2 + b^2) - bi/(a^2 + b^2) */
		  if(!tmp1 && !tmp2) {
			fdpfilterdata[index].re = 0;
			fdpfilterdata[index].im = 0;
		  }
		  else {
			fdpfilterdata[index].re = tmp1 / ( pow(tmp1,2)+pow(tmp2,2) ) * scale;
			fdpfilterdata[index].im = -tmp2 / ( pow(tmp1,2)+pow(tmp2,2) ) * scale;
		  }
		  tmp1 = sqrt(pow(fdpfilterdata[index].im,2)+pow(fdpfilterdata[index].re,2));
		  if(scale == 1.0 &&  tmp1 < 1.0) {
			fprintf(stderr, "BOOP BOOP ERROR!  fdpfilter mag == %f at value z,y,x=%d,%d,%d.\n",tmp1, k,j,i);
			fprintf(stderr, "psfindex1=%d, psfindex2=%d, , depth = %f. idepth = %d, frac=%f.\n", 
					psfindex1, psfindex2, depth, idepth, frac);
			fprintf(stderr, "fdr[psfindex1]= %f + %f i, fdr[psfindex2] = %f + %f i\n", 
					psffftdata[psfindex1].re, psffftdata[psfindex1].im,
					psffftdata[psfindex2].re, psffftdata[psfindex2].im);
			fprintf(stderr, "Slope is %f, scale is %f\n", slope, scale);
			fdpfilterdata[index].re = fdpfilterdata[index].re * 1./tmp1;
			fdpfilterdata[index].im = fdpfilterdata[index].im * 1./tmp1;
		  }
		}
		else {
		  fdpfilterdata[index].re = tmp1 * scale;
		  fdpfilterdata[index].im = tmp2 * scale;
		}
	  }
	}
  }
  
		/*
		if(flags == INVERSE_FILTER) {




		  //  1/(a+bi) = a/(a^2 +b^2) - bi/(a^2+b^2) 
		  // Calculate demoninators in order to save time.  Ha! 
		  denom1 = (psffftdata[psfindex1].re*psffftdata[psfindex1].re +\
					psffftdata[psfindex1].im*psffftdata[psfindex1].im);
		  denom2 = (psffftdata[psfindex2].re*psffftdata[psfindex2].re +\
					psffftdata[psfindex2].im*psffftdata[psfindex2].im);
		  // *************
		  // * Real first 
		  // *************
		  tmp1 = psffftdata[psfindex1].re / denom1;
		  tmp2 = psffftdata[psfindex2].re / denom2;
		  // Want to avoid nan 
		  if(!psffftdata[psfindex1].re && !psffftdata[psfindex1].im)
			tmp1 = 0;
		  if(!psffftdata[psfindex2].re && !psffftdata[psfindex2].im)
			tmp2 = 0;
		  // Now interpolate for final result 
		  fdpfilterdata[index].re = ((1-frac)*tmp1 + frac*tmp2)*scale;
		  // ****************
		    * Now imaginary
		    ****************
		  tmp1 = -psffftdata[psfindex1].im / denom1;
		  tmp2 = -psffftdata[psfindex2].im / denom2;
		  // Want to avoid nan 
		  if(!psffftdata[psfindex1].re && !psffftdata[psfindex1].im)
			tmp1 = 0;
		  if(!psffftdata[psfindex2].re && !psffftdata[psfindex2].im)
			tmp2 = 0;
		  // Now do interpolation for final result 
		  fdpfilterdata[index].im = ((1-frac)*tmp1 + frac*tmp2)*scale;
		}
		else {
		  //		  fprintf(stdout, "Values are %f %f\n", psffftdata[psfindex1].re, psffftdata[psfindex1].im);
		  fdpfilterdata[index].re = psffftdata[psfindex1].re*scale;
		  fdpfilterdata[index].im = psffftdata[psfindex1].im*scale;
		}
	  } 
	  }
	*/ 
	//fprintf(stdout, "First value is %f + %f.i\n", fdpfilterdata[index].re, fdpfilterdata[index].im);
  if(verbose) fprintf(stdout,"Completed building the FDP filter.\n");
  
  /* Finally, we're successful */
  if(verbose) fprintf(stdout, "*********** Returning from build_fdp **************\n");
  
  return 1;
}

/***************************************************************
 ***************************************************************
 ******************* DATA MANIPULATION *************************
 ***************************************************************
 ***************************************************************/



int multiply_complex_by_complex(long nelems, fdr_complex *data1, fdr_complex *data2)
{
  long i;
  float x, y;

  //#pragma omp parallel for
  for(i=0;i<nelems;i++) {
	/* (a+bi)*(c+di) = (ac-bd)+(ad+bc)i */
	if(i==151901747) {
	  fprintf(stdout, "At index %d, data1 = %f + %f i , data2 = %f + %f i\n", i, data1[i].re, data1[i].im, data2[i].re, data2[i].im);
	  fprintf(stdout, "we are expecting %f + %f i\n", data1[i].re*data2[i].re - data1[i].im*data2[i].im,data1[i].re*data2[i].im+data1[i].im*data2[i].re); 
	}
	
	x = data1[i].re*data2[i].re - data1[i].im*data2[i].im;
	y = data1[i].re*data2[i].im + data1[i].im*data2[i].re;
	data1[i].re = x;
	data1[i].im = y;
	if(i==151901747) 
	  fprintf(stdout, "At index %d, result is data1 = %f + %f i\n", i, data1[i].re, data1[i].im);

  }
  return 1;
}

int multiply_complex_by_real(long nelems, fdr_complex *data1, float *data2)
{
  long i;

  /* k*(a+bi) = ka + kbi */
  //#pragma omp parallel for
  for(i=0;i<nelems;i++) {
	data1[i].re = data1[i].re*data2[i];
	data1[i].im = data1[i].im*data2[i];
  }
  return 1;
}

int multiply_complex_by_gaussian(long n3, long n2, long n1, fdr_complex *data, double weight)
{
  long k,j,i,index;
  float *gauss3, *gauss2, *gauss1;
  float n32, n22, n12;
  
  gauss3 = (float *)calloc(n3,sizeof(float));
  gauss2 = (float *)calloc(n2,sizeof(float));
  gauss1 = (float *)calloc(n1,sizeof(float));
  n32 = (n3-1)/2.;
  n22 = (n2-1)/2.;
  n12 = (n1-1)/1.;
  
  if(gauss3 == NULL || gauss2 == NULL || gauss1 == NULL) {
	fprintf(stderr, "Could not allocate memory for gaussians!\n");
	free(gauss3); free(gauss2); free(gauss1);
	return 0;
  }

  /* First let's prepare the gaussian in each coordinate to minimize calcs  */ 
  for(k=0;k<n3;k++){ 
	gauss3[k] = exp(-(k-n32)*(k-n32)/(2*weight*n32)*(2*weight*n32));
  }
  for(j=0;j<n2;j++) {
	gauss2[j] = exp(-(j-n22)*(j-n22)/(2*weight*n22)*(2*weight*n22));
  }
  for(i=0;i<n1;i++) {
	gauss1[i] = exp(-(i-n12)*(i-n12)/(2*weight*n12)*(2*weight*n12));
  }
  
#pragma omp parallel for private(index,i,j)
  for(k=0;k<n3;k++) {
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		data[index].re = data[index].re * gauss3[k]*gauss2[j]*gauss1[i];
		data[index].im = data[index].im * gauss3[k]*gauss2[j]*gauss1[i];
	  }
	}
  }

  return 1;

}

int multiply_real_by_gaussian(long n3, long n2, long n1, float *data, double wght)
{
  long k,j,i,index;
  float *gauss3, *gauss2, *gauss1;
  float n32, n22, n12;
  
  gauss3 = (float *)calloc(n3,sizeof(float));
  gauss2 = (float *)calloc(n2,sizeof(float));
  gauss1 = (float *)calloc(n1,sizeof(float));
  n32 = (n3-1)/2.;
  n22 = (n2-1)/2.;
  n12 = (n1-1)/2.;

  if(gauss3 == NULL || gauss2 == NULL || gauss1 == NULL) {
	fprintf(stderr, "Could not allocate memory for gaussians!\n");
	free(gauss3); free(gauss2); free(gauss1);
	return 0;
  }

  /* First let's prepare the gaussian in each coordinate to minimize calcs  */ 
  for(k=0;k<n3;k++) {
	gauss3[k] = exp(- ((k-n32)*(k-n32)) / ((2*wght*n32)*(2*wght*n32)) );
  }
  for(j=0;j<n2;j++) {
	gauss2[j] = exp(- ((j-n22)*(j-n22)) / ((2*wght*n22)*(2*wght*n22)) );
  }
  for(i=0;i<n1;i++) {
	gauss1[i] = exp(- ((i-n12)*(i-n12)) / ((2*wght*n12)*(2*wght*n12)) );
  }
  
#pragma omp parallel for private(j,i,index)
  for(k=0;k<n3;k++) {
	for(j=0;j<n2;j++) {
	  for(i=0;i<n1;i++) {
		index = k*n2*n1 + j*n1 + i;
		//fprintf(stdout, "Multiplying %f by 1) %f 2) %f 3) %f\n", data[index],gauss3[k],gauss2[k],gauss1[j]);
		data[index] = data[index] * gauss3[k]*gauss2[j]*gauss1[i];
	  }
	}
  }

  return 1;

}
