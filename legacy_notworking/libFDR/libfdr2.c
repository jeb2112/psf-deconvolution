#include "fdr2.h"
#define NUMSLICES 10

int verbose = 2;



int shift_float_to_complex(long elems, fdr_complex *data) {
  float *tmp;
  long i;

  tmp = (float *)data;
  tmp += elems-1;
  // Start at the end of the data and transpose
  for(i=elems-1;i>=0;i--) 
  {
		// fprintf(stdout,"Index is %d\n",i);
		data[i]= (float complex) *tmp;
		// fprintf(stdout,"pointer value is %f, data is %f %f\n", *tmp, creal(data[i]), cimag(data[i]));
		tmp = tmp - 1;
  }

  return 1;
}


/***************************************************************
 ***************************************************************
 ****************** MINC FILE FUNCTIONS ************************
 ***************************************************************
 ***************************************************************/


int get_minc_dimensions_from_handle(mihandle_t filein, long *nz, long *ny, long *nx)
{
  midimhandle_t in_dimensions[3];
  unsigned int in_dimsizes[3];
  int result;

  /* Get the dimension handles */
  result = miget_volume_dimensions (filein, MI_DIMCLASS_SPATIAL, 
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


int open_minc_file_read(char *filename, mihandle_t *handle, int flags)
{
  int result;

  /* Get volume handle from minc file */
  result = miopen_volume(filename,MI2_OPEN_READ, handle);
  if (result == MI_ERROR) 
  {
		fprintf(stderr, "Error opening the input file %s.\n", filename);
		return 0;
  } 
  else if(verbose>1) fprintf(stdout,"Opened file\n");
	
	return 1;
	
}

int open_minc_file_write(char *filename, mihandle_t *handle, long n3, long n2, long n1, int flags)
{
 	mihandle_t out_volume;
  midimhandle_t out_dimensions[3];
  unsigned int out_dimsizes[3];
  unsigned long out_count[3];
  unsigned long start[3];
  long nelems;

  float tmp, tmp2;
  int result;
  
  if(verbose>1) 
  	fprintf(stdout, "Entering write minc file\n");
  
  start[0] = start[1] = start[2] = 0;
  out_dimsizes[0] = n3;
  out_dimsizes[1] = n2;
  out_dimsizes[2] = n1;
  nelems = n3*n2*n1;

  /* Create the three dimensions */
  result = micreate_dimension("xspace",MI_DIMCLASS_SPATIAL,MI_DIMATTR_REGULARLY_SAMPLED,
							  out_dimsizes[2],&out_dimensions[2]);
  if(result == MI_ERROR) 
  { 
		fprintf(stderr, "Error creating output dimensions\n"); 
		return 0;
  }  
  result = micreate_dimension("yspace",MI_DIMCLASS_SPATIAL,MI_DIMATTR_REGULARLY_SAMPLED,
							  out_dimsizes[1],&out_dimensions[1]);
  if(result == MI_ERROR) 
  { 
		fprintf(stderr, "Error creating output dimensions\n"); 
		return 0;
  }  
  result = micreate_dimension("zspace",MI_DIMCLASS_SPATIAL,MI_DIMATTR_REGULARLY_SAMPLED,
							  out_dimsizes[0],&out_dimensions[0]);
  if(result == MI_ERROR) 
  { 
		fprintf(stderr, "Error creating output dimensions\n"); 
		return 0;
  }  
  else if(verbose>1) fprintf(stdout, "Created dimensions for output file.\n");

	
  /* Create the file on the disk and get volume handle*/
  if(flags == COMPLEX)
		result = micreate_volume(filename, 3, out_dimensions,MI_TYPE_FCOMPLEX, MI_CLASS_COMPLEX,NULL,handle);
	else
		result = micreate_volume(filename, 3, out_dimensions, MI_TYPE_FLOAT, MI_CLASS_REAL,NULL,handle);
  if(result == MI_ERROR) 
  {
		fprintf(stderr, "Error creating output volume\n");
		return 0;
  }
  else if(verbose>1) fprintf(stdout, "Created volume for output.\n");


	/* Create the image associated with the volume */
	result = micreate_volume_image(*handle);
	if(result == MI_ERROR)
	{
		fprintf(stderr, "Error creating volume image.\n");
		return 0;
	}

  return 1;
	
}

/***************************************************************
 ***************************************************************
 ******************* DATA COMPARE FUNCTIONS ********************
 ***************************************************************
 ***************************************************************/


int slicewise_compare_complex_data(mihandle_t file1, mihandle_t file2) {

  long nz,ny,nx,mz,my,mx,i,k;
  unsigned long start[3], count[3];
  int result;
  fdr_complex *data1, *data2;
  
  result = get_minc_dimensions_from_handle(file1, &mz, &my, &mx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim1\n");
  
  result = get_minc_dimensions_from_handle(file2, &nz, &ny, &nx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim2\n");
  
  if(mz != nz || my != ny || mx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in1: [%d,%d,%d] in2: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  
  data1 = (fdr_complex *)fftwf_malloc(ny*nx*sizeof(fdr_complex));
  data2 = (fdr_complex *)fftwf_malloc(ny*nx*sizeof(fdr_complex));
  
  fprintf(stdout, "Beginning comparison...\n");
  for(k=0;k<nz;k++) {
    start[0] = (unsigned long)k; start[1] = start[2] = 0;
    count[0] = 1;
    count[1] = (unsigned long)ny;
    count[2] = (unsigned long)nx;
    
    result = miget_voxel_value_hyperslab(file1,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) data1);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading data1.\n");
      return 0;
    }
    
    result = miget_voxel_value_hyperslab(file2,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) data2);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading data2.\n");
      return 0;
    }
    for(i=0;i<ny*nx;i++) {
      if(creal(data1[i]) != creal(data2[i]))
        fprintf(stdout, "At k=%d,i=%d, REAL data1 is different than data2 (%f vs %f)\n", k,i,creal(data1[i]), creal(data2[i]));
      if(cimag(data1[i]) != cimag(data2[i]))
        fprintf(stdout, "At k=%d,i=%d, REAL data1 is larger than data2 (%f vs %f)\n", k,i,cimag(data1[i]),cimag(data2[i]));
    }
  }
  
  fftwf_free(data1); fftwf_free(data2);
  
  fprintf(stdout, "Ending comparison.  If you don't see any errors, you're good!\n");

  return 1;

}

//int compare_real_data(long n3,long n2,long n1, float *data1, float *data2) {
//
//  long j,k,i, index;
//  
//  fprintf(stdout, "Beginning comparison...\n");
//  for(k=0;k<n3;k++) {
//  for(j=0;j<n2;j++) {
//    for(i=0;i<n1;i++) {
//    index = k*n2*n1 + j*n1 + i;
//    if(data1[index] > data2[index]) {
//      fprintf(stdout, "At k=%d,j=%d,i=%d, data1 is larger than data2 (%f vs %f)\n", k,j,i,data1[index], data2[index]);
//    }
//    }
//  }
//  }
//  fprintf(stdout, "Ending comparison.  If you don't see any errors, you're good!\n");
//
//  return 1;
//
//}



/***************************************************************
 ***************************************************************
 ************** SLICEWISE MINC MANIPULATION ********************
 ***************************************************************
 ***************************************************************/
 
 int slicewise_rfftshift(mihandle_t filein, mihandle_t fileout) {
  
  float tmp;
  long i,j,k, index1, index2, Z, Y, X, Y2, Z2, X2, mz, my, mx, nz1, ny1, nx1;
  float *temp2d_q1 = NULL; float *temp2d_q2 = NULL;
  unsigned long start_q1[3], start_q2[3], count[3];
  
  int result;
  
  fprintf(stdout, "Entering fftshift\n");
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim1\n");
  
  result = get_minc_dimensions_from_handle(fileout, &nz1, &ny1, &nx1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim2\n");
  
  if(mz != nz1 || my != ny1 || mx != nx1)
  {
    fprintf(stderr, "Dimensions do not match (in1: [%d,%d,%d] in2: [%d,%d,%d])\n", mz,my,mx, nz1,ny1,nx1);
    return 0;
  }
  
  Z = nz1; Y = ny1; X = nx1;
  Z2 = Z/2;
  Y2 = Y/2;
  X2 = X/2;
  temp2d_q1 = (float *)fftwf_malloc(Y*X*sizeof(float));
  temp2d_q2 = (float *)fftwf_malloc(Y*X*sizeof(float));
  for(k=0;k<Z2; k++) {
    
    fprintf(stdout, "FFTSHIFT Z=%d\n", k);
    
    /* Load the data in 2 slices */
    start_q1[0] = (unsigned long)k; start_q1[1] = start_q1[2] = 0;
    start_q2[0] = (unsigned long)k+Z2; start_q2[1] = start_q2[2] = 0;
    count[0] = 1; count[1] = (unsigned long)ny1; count[2] = (unsigned long)nx1;
    
    //fprintf(stdout, "Getting from %d,%d,%d and %d,%d,%d for %d,%d,%d\n", 
    //        start_q1[0], start_q1[1], start_q1[2], start_q2[0], start_q2[1],
    //       start_q2[2], count[0], count[1], count[2]);
    result = miget_voxel_value_hyperslab(filein,MI_TYPE_FLOAT,start_q1,count, (float *) temp2d_q1);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading slice1.\n");
      return 0;
    }
    
    result = miget_voxel_value_hyperslab(filein,MI_TYPE_FLOAT,start_q2,count, (float *) temp2d_q2);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading slice 2.\n");
      return 0;
    }
    
    /* Data from two slices is now loaded, so do the swap */
    for(j=0;j<Y2;j++) {
      for(i=0;i<X2;i++) {
      // Swap quadrant 1 and 5
        index1 = j*X + i;
        index2 = (j+Y2)*X + i+X2;
        tmp = temp2d_q1[index1];
        temp2d_q1[index1] = temp2d_q2[index2];
        temp2d_q2[index2] = tmp;
      // Swap quadrant 2 and 6
        index1 = j*X + i+X2;
        index2 = (j+Y2)*X + i;
        tmp = temp2d_q1[index1];
        temp2d_q1[index1] = temp2d_q2[index2];
        temp2d_q2[index2] = tmp;
      // Swap quadrant 3 and 7
        index1 = (j+Y2)*X + i;
        index2 = j*X + i+X2;
        tmp = temp2d_q1[index1];
        temp2d_q1[index1] = temp2d_q2[index2];
        temp2d_q2[index2] = tmp;
      // Swap quadrant 4 and 8
        index1 = (j+Y2)*X + i+X2;
        index2 = j*X + i;
        tmp = temp2d_q1[index1];
        temp2d_q1[index1] = temp2d_q2[index2];
        temp2d_q2[index2] = tmp;
      } // end xloop
    } // end yloop
    
    /* Write out the data */
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FLOAT,start_q1,count, (float *) temp2d_q1);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error writing slice 1.\n");
      return 0;
    }
    
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FLOAT,start_q2,count, (float *) temp2d_q2);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error writing slice 2.\n");
      return 0;
    }
    
  } // end zloop

  fftwf_free(temp2d_q1); fftwf_free(temp2d_q2);
  
  return 1;
}
 
int slicewise_cfftshift(mihandle_t filein, mihandle_t fileout) {
  
  fdr_complex tmp;
  long i,j,k, index1, index2, Z, Y, X, Y2, Z2, X2, mz, my, mx, nz, ny, nx;
  fdr_complex *temp2d_q1 = NULL; fdr_complex *temp2d_q2 = NULL;
  unsigned long start_q1[3], start_q2[3], count[3];
  
  int result;
  
  fprintf(stdout, "Entering fftshift\n");
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim1\n");
  
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim2\n");
  
  if(mz != nz || my != ny || mx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in1: [%d,%d,%d] in2: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  
  Z = nz; Y = ny; X = nx;
  Z2 = Z/2;
  Y2 = Y/2;
  X2 = X/2;
  temp2d_q1 = (fdr_complex *)fftwf_malloc(Y*X*sizeof(fdr_complex));
  temp2d_q2 = (fdr_complex *)fftwf_malloc(Y*X*sizeof(fdr_complex));
  for(k=0;k<Z2; k++) {
    
    fprintf(stdout, "FFTSHIFT Z=%d\n", k);
    
    /* Load the data in 2 slices */
    start_q1[0] = (unsigned long)k; start_q1[1] = start_q1[2] = 0;
    start_q2[0] = (unsigned long)k+Z2; start_q2[1] = start_q2[2] = 0;
    count[0] = 1; count[1] = (unsigned long)ny; count[2] = (unsigned long)nx;
    
    fprintf(stdout, "Getting from %d,%d,%d and %d,%d,%d for %d,%d,%d\n", 
            start_q1[0], start_q1[1], start_q1[2], start_q2[0], start_q2[1],
           start_q2[2], count[0], count[1], count[2]);
    result = miget_voxel_value_hyperslab(filein,MI_TYPE_FCOMPLEX,start_q1,count, (fdr_complex *) temp2d_q1);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading slice1.\n");
      return 0;
    }
    
    result = miget_voxel_value_hyperslab(filein,MI_TYPE_FCOMPLEX,start_q2,count, (fdr_complex *) temp2d_q2);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading slice 2.\n");
      return 0;
    }
    
    /* Data from two slices is now loaded, so do the swap */
    for(j=0;j<Y2;j++) {
      for(i=0;i<X2;i++) {
      // Swap quadrant 1 and 5
        index1 = j*X + i;
        index2 = (j+Y2)*X + i+X2;
        tmp = temp2d_q1[index1];
        temp2d_q1[index1] = temp2d_q2[index2];
        temp2d_q2[index2] = tmp;
      // Swap quadrant 2 and 6
        index1 = j*X + i+X2;
        index2 = (j+Y2)*X + i;
        tmp = temp2d_q1[index1];
        temp2d_q1[index1] = temp2d_q2[index2];
        temp2d_q2[index2] = tmp;
      // Swap quadrant 3 and 7
        index1 = (j+Y2)*X + i;
        index2 = j*X + i+X2;
        tmp = temp2d_q1[index1];
        temp2d_q1[index1] = temp2d_q2[index2];
        temp2d_q2[index2] = tmp;
      // Swap quadrant 4 and 8
        index1 = (j+Y2)*X + i+X2;
        index2 = j*X + i;
        tmp = temp2d_q1[index1];
        temp2d_q1[index1] = temp2d_q2[index2];
        temp2d_q2[index2] = tmp;
      } // end xloop
    } // end yloop
    
    /* Write out the data */
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start_q1,count, (mifcomplex_t *) temp2d_q1);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error writing slice 1.\n");
      return 0;
    }
    
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start_q2,count, (mifcomplex_t *) temp2d_q2);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error writing slice 2.\n");
      return 0;
    }
    
  } // end zloop

  fftwf_free(temp2d_q1); fftwf_free(temp2d_q2);
  
  return 1;
}

int slicewise_multiply_complex_by_real(mihandle_t file1, mihandle_t file2, mihandle_t fileout, int flags)
{
  long mz,my,mx;
  long nz,ny,nx;
  long pz,py,px;
  long z, slice, i;
  unsigned long start[3], count[3];
  
  fdr_complex *data1_3d = NULL;
  float *data2_3d = NULL;
  fdr_complex *data3_3d = NULL;
  
  float x, y;

  int result;
  
  result = get_minc_dimensions_from_handle(file1, &mz, &my, &mx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim1\n");
  
  result = get_minc_dimensions_from_handle(file2, &nz, &ny, &nx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim2\n");
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in1: [%d,%d,%d] in2: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  
  result = get_minc_dimensions_from_handle(fileout, &pz, &py, &px);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dimout\n");
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in1: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }  
  
  data1_3d = (fdr_complex *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(fdr_complex));
  data2_3d = (float *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(float));
  data3_3d = (fdr_complex *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(fdr_complex));
  if(verbose>1) fprintf(stdout, "Allocated data\n");
  
  for(z=0;z<nz;z+=NUMSLICES)
  {   
    /* 
     * Retrieve the data from the file 
     */
    start[0] = (unsigned long)z; start[1] = start[2] = 0;
    count[0] = NUMSLICES;
    count[1] = (unsigned long)my;
    count[2] = (unsigned long)mx;
    
    fprintf(stdout, "Collecting data from %d %d %d -- %d %d %d\n", start[0], start[1], start[2],
            count[0]+start[0], count[1], count[2]);
    result = miget_voxel_value_hyperslab(file1,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) data1_3d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading 2D data file1.\n");
      return 0;
    }
    
    result = miget_voxel_value_hyperslab(file2,MI_TYPE_FLOAT,start,count, (float *) data2_3d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading 2D data file2.\n");
      return 0;
    }
        
    for(i=0;i<nx*ny*NUMSLICES;i++)
    {
      // fprintf(stdout, "Inputs are %2.3f + %2.3fi * %2.3f", creal(data1_3d[i]), cimag(data1_3d[i]), data2_3d[i]);
      data3_3d[i] = ( (float)creal(data1_3d[i]) * (float)data2_3d[i] ) + ( (float)cimag(data1_3d[i]) * (float)data2_3d[i] ) * I;
      // fprintf(stdout, "And the output is %2.3f + %2.3fi\n", creal(data3_3d[i]) cimag(data3_3d[i]));
    }
    
    /* Write that out to the file */
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)data3_3d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error in setting values.\n");
      return 0;
    }
  }
  fftwf_free(data1_3d); fftwf_free(data2_3d); fftwf_free(data3_3d);
  
  return 1;
}

int slicewise_multiply_complex_by_complex(mihandle_t file1, mihandle_t file2, mihandle_t fileout, int flags)
{
  
  long mz,my,mx;
  long nz,ny,nx;
  long pz,py,px;
  long z, slice, i;
  unsigned long start[3], count[3];
  
  fdr_complex *data1_3d = NULL;
  fdr_complex *data2_3d = NULL;
  fdr_complex *data3_3d = NULL;
  
  float x, y;

  int result;
  
  result = get_minc_dimensions_from_handle(file1, &mz, &my, &mx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim1\n");
  
  result = get_minc_dimensions_from_handle(file2, &nz, &ny, &nx);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dim2\n");
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in1: [%d,%d,%d] in2: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  
  result = get_minc_dimensions_from_handle(fileout, &pz, &py, &px);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "Obtained dimout\n");
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in1: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }  
  
  data1_3d = (fdr_complex *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(fdr_complex));
  data2_3d = (fdr_complex *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(fdr_complex));
  data3_3d = (fdr_complex *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(fdr_complex));
  if(verbose>1) fprintf(stdout, "Allocated data\n");
  
  for(z=0;z<nz;z+=NUMSLICES)
  {   
    /* 
     * Retrieve the data from the file 
     */
    start[0] = (unsigned long)z; start[1] = start[2] = 0;
    count[0] = NUMSLICES;
    count[1] = (unsigned long)my;
    count[2] = (unsigned long)mx;
    
    fprintf(stdout, "Collecting data from %d %d %d -- %d %d %d\n", start[0], start[1], start[2],
            count[0]+start[0], count[1], count[2]);
    result = miget_voxel_value_hyperslab(file1,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) data1_3d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading 2D data file1.\n");
      return 0;
    }
    
    result = miget_voxel_value_hyperslab(file2,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) data2_3d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error loading 2D data file2.\n");
      return 0;
    }
        
    for(i=0;i<nx*ny*NUMSLICES;i++)
    {
      data3_3d[i] = (fdr_complex)(data1_3d[i] * data2_3d[i]);
    }
    
    /* Write that out to the file */
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)data3_3d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error in setting values.\n");
      return 0;
    }
  }
  fftwf_free(data1_3d); fftwf_free(data2_3d); fftwf_free(data3_3d);
  
  return 1;
}


/***************************************************************
 ***************************************************************
 ****************** FULL DATA FUNCTIONS ************************
 ***************************************************************
 ***************************************************************/
 
int full_magnitude(mihandle_t filein, mihandle_t fileout, int flags)
{
  long mz,my,mx;
  long nz,ny,nx;
  long slice, x,y,z, i;
  fdr_complex *datain = NULL;
  fdr_complex *dataout = NULL;
  fdr_complex tmp;
  float tmp2;
  unsigned long start[3], count[3],nelems;
  float scale;
  double dmin,dmax;

  int result;
  
  
  dmin = 1e100; dmax = 0;
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  
  // Dimensions are the same size, we can feel free to do the transforms!

  datain = (fdr_complex *)fftwf_malloc(nz*ny*nx*sizeof(fdr_complex));
  dataout = (fdr_complex *)fftwf_malloc(nz*ny*nx*sizeof(float));
  
  start[0] = (unsigned long)0; start[1] = start[2] = 0;
  count[0] = (unsigned long)mz;
  count[1] = (unsigned long)my;
  count[2] = (unsigned long)mx;
  result = miget_voxel_value_hyperslab(filein,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) datain);
  if(result == MI_ERROR) {
    fprintf(stderr, "Error getting data.\n");
    return 0;
  }
  else if(verbose>1) fprintf(stdout,"Retrieved data\n");
  
  
  /*  Calculate the fabs  */
  for(i=0;i<nx*ny*nz;i++) {
    tmp2 = (float)fabs(datain[i]);
    if(tmp2 < dmin)
      dmin = tmp2;
    if(tmp2 > dmax)
      dmax = tmp2;
    if(isinf(tmp2)) 
      fprintf(stdout, "Infinity calculated at index %d, from %f %f\n", i, creal(datain[i]), creal(datain[i]));
  }
  
  /* Write that out to the file */
  result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FLOAT,start,count,(float *)dataout);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Error in setting values.\n");
    return 0;
  }

  fftwf_free(datain); fftwf_free(dataout); 
  
  if(verbose)
    fprintf(stdout, "Completed full magnitude\n");
    
  return 1;
}

int series_1dfft(mihandle_t filein, mihandle_t fileout, int flags)
{
  long mz,my,mx;
  long nz,ny,nx;
  long slice, x,y,z, i;
  fdr_complex *datain = NULL;
  fdr_complex *dataout = NULL;
  fftwf_plan plan1d;
  unsigned long start[3], count[3],nelems;
  float scale;
  
  int result;
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  if(my != 1)
  {
    fprintf(stderr, "Yspace must have dimensions of 1!!\n");
    return 0;
  }
  
  // Dimensions are the same size, we can feel free to do the transforms!

  datain = (fdr_complex *)fftwf_malloc(nx*sizeof(fdr_complex));
  dataout = (fdr_complex *)fftwf_malloc(nx*sizeof(fdr_complex));
  
  plan1d = fftwf_plan_dft_1d(mx, datain, dataout, FFTW_FORWARD, FFTW_MEASURE);
  
  for(z=0;z<nz;z++)
  {
    start[0] = (unsigned long)z; start[1] = start[2] = 0;
    count[0] = (unsigned long)1;
    count[1] = (unsigned long)my;
    count[2] = (unsigned long)mx;
    result = miget_voxel_value_hyperslab(filein,MI_TYPE_FLOAT,start,count, (float *) datain);
    if(result == MI_ERROR) {
      fprintf(stderr, "Error getting data.\n");
      return 0;
    }
    else if(verbose>1) fprintf(stdout,"Retrieved data\n");
    
    shift_float_to_complex(nx, (float *) datain);
    if(verbose>1) fprintf(stdout, "Rearranged data from float to complex\n");
    
    fftwf_execute(plan1d);
    if(verbose>1) fprintf(stdout, "Executed 1D FFT at z=%d\n",z);
    
    /* Now scale the data since FFTW does not*/
    scale = 1.0/sqrt(1.0*nx);
    for(i=0;i<nx;i++)
      dataout[i] *= scale;
    
    /* Write that out to the file */
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)dataout);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error in setting values.\n");
      return 0;
    }
  }
  
  fftwf_destroy_plan(plan1d);   
  fftwf_free(datain); fftwf_free(dataout); 
  
  fprintf(stdout, "Right, so we've finished the z stack, let's do the x stack.\n");
  
  /* Now in the other direction */
  datain = (fdr_complex *)fftwf_malloc(nz*sizeof(fdr_complex));
  dataout = (fdr_complex *)fftwf_malloc(nz*sizeof(fdr_complex));
  
  plan1d = fftwf_plan_dft_1d(mz, datain, dataout, FFTW_FORWARD, FFTW_MEASURE);
  
  for(x=0;x<nx;x++)
  {
    start[0] = 0; start[1] = 0; start[2] = (unsigned long)x;
    count[0] = (unsigned long)mz;
    count[1] = (unsigned long)my;
    count[2] = (unsigned long)1;
    result = miget_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) datain);
    if(result == MI_ERROR) {
      fprintf(stderr, "Error getting data.\n");
      return 0;
    }
    else if(verbose>1) fprintf(stdout,"Retrieved data\n");

    fftwf_execute(plan1d);
    if(verbose>1) fprintf(stdout, "Executed 1D FFT at x=%d\n",x);
    
    /* Now scale the data since FFTW does not*/
    scale = 1.0/sqrt(1.0*nz);
    for(i=0;i<nz;i++)
      dataout[i] *= scale;
    
    /* Write that out to the file */
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)dataout);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error in setting values.\n");
      return 0;
    }
  }
  fftwf_destroy_plan(plan1d);   
  fftwf_free(datain); fftwf_free(dataout); 
   
  if(verbose)
    fprintf(stdout, "Completed 2nd series of 1D FFTs\n");
    
  return 1;
}

int full_2dfft(mihandle_t filein, mihandle_t fileout, int flags)
{
  long mz,my,mx;
  long nz,ny,nx;
  long slice, x,y,z, i;
  fdr_complex *datain = NULL;
  fdr_complex *dataout = NULL;
  fftwf_plan plan2d;
  unsigned long start[3], count[3],nelems;
  float scale;
  
  int result;
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  if(my != 1)
  {
    fprintf(stderr, "Yspace must have dimensions of 1!!\n");
    return 0;
  }
  
  // Dimensions are the same size, we can feel free to do the transforms!

  datain = (fdr_complex *)fftwf_malloc(nz*nx*sizeof(fdr_complex));
  dataout = (fdr_complex *)fftwf_malloc(nz*nx*sizeof(fdr_complex));
  
  plan2d = fftwf_plan_dft_2d(mz,mx, datain, dataout, FFTW_FORWARD, FFTW_MEASURE);
  
  start[0] = (unsigned long)0; start[1] = start[2] = 0;
  count[0] = (unsigned long)mz;
  count[1] = (unsigned long)my;
  count[2] = (unsigned long)mx;
  result = miget_voxel_value_hyperslab(filein,MI_TYPE_FLOAT,start,count, (float *) datain);
  if(result == MI_ERROR) {
    fprintf(stderr, "Error getting data.\n");
    return 0;
  }
  else if(verbose>1) fprintf(stdout,"Retrieved data\n");
  
  shift_float_to_complex(nz*ny*nx, (float *) datain);
  if(verbose>1) fprintf(stdout, "Rearranged data from float to complex\n");
  
  fftwf_execute(plan2d);
  if(verbose>1) fprintf(stdout, "Executed 2D FFT\n");
  
  /* Now scale the data since FFTW does not*/
  scale = 1.0/sqrt(1.0*nz*nx*ny);
  for(i=0;i<nz*ny*nx;i++)
    dataout[i] *= scale;
  
  /* Write that out to the file */
  result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)dataout);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Error in setting values.\n");
    return 0;
  }
  
  fftwf_destroy_plan(plan2d);   
  fftwf_free(datain); fftwf_free(dataout); 
  
  if(verbose)
    fprintf(stdout, "Completed full 2d FFT\n");
    
  return 1;
}

int full_2difft(mihandle_t filein, mihandle_t fileout, int flags)
{
  long mz,my,mx;
  long nz,ny,nx;
  long slice, x,y,z, i;
  fdr_complex *datain = NULL;
  fdr_complex *dataout = NULL;
  fftwf_plan plan2d;
  unsigned long start[3], count[3],nelems;
  float scale;
  
  int result;
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  if(my != 1)
  {
    fprintf(stderr, "Yspace must have dimensions of 1!!\n");
    return 0;
  }
  
  // Dimensions are the same size, we can feel free to do the transforms!

  datain = (fdr_complex *)fftwf_malloc(nz*ny*nx*sizeof(fdr_complex));
  dataout = (fdr_complex *)fftwf_malloc(nz*ny*nx*sizeof(fdr_complex));
  
  plan2d = fftwf_plan_dft_3d(mz,my,mx, datain, dataout, FFTW_BACKWARD, FFTW_MEASURE);
  
  start[0] = (unsigned long)0; start[1] = start[2] = 0;
  count[0] = (unsigned long)mz;
  count[1] = (unsigned long)my;
  count[2] = (unsigned long)mx;
  result = miget_voxel_value_hyperslab(filein,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) datain);
  if(result == MI_ERROR) {
    fprintf(stderr, "Error getting data.\n");
    return 0;
  }
  else if(verbose>1) fprintf(stdout,"Retrieved data\n");
  
  fftwf_execute(plan2d);
  if(verbose>1) fprintf(stdout, "Executed 2D FFT\n");
  
  /* Now scale the data since FFTW does not*/
  scale = 1.0/sqrt(1.0*nz*nx*ny);
  for(i=0;i<nz*ny*nx;i++)
    dataout[i] *= scale;
  
  /* Write that out to the file */
  result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)dataout);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Error in setting values.\n");
    return 0;
  }
  
  fftwf_destroy_plan(plan2d);
  fftwf_free(datain); fftwf_free(dataout); 
  
  if(verbose)
    fprintf(stdout, "Completed full 2d FFT\n");
    
  return 1;
}

int full_3dfft(mihandle_t filein, mihandle_t fileout, int flags)
{
  long mz,my,mx;
  long nz,ny,nx;
  long slice, x,y,z, i;
  fdr_complex *datain = NULL;
  fdr_complex *dataout = NULL;
  fftwf_plan plan3d;
  unsigned long start[3], count[3],nelems;
  float scale;
  
  int result;
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  
  // Dimensions are the same size, we can feel free to do the transforms!

  datain = (fdr_complex *)fftwf_malloc(nz*ny*nx*sizeof(fdr_complex));
  dataout = (fdr_complex *)fftwf_malloc(nz*ny*nx*sizeof(fdr_complex));
  
  plan3d = fftwf_plan_dft_3d(mz,my,mx, datain, dataout, FFTW_FORWARD, FFTW_MEASURE);
  
  start[0] = (unsigned long)0; start[1] = start[2] = 0;
  count[0] = (unsigned long)mz;
  count[1] = (unsigned long)my;
  count[2] = (unsigned long)mx;
  result = miget_voxel_value_hyperslab(filein,MI_TYPE_FLOAT,start,count, (float *) datain);
  if(result == MI_ERROR) {
    fprintf(stderr, "Error getting data.\n");
    return 0;
  }
  else if(verbose>1) fprintf(stdout,"Retrieved temp3d data\n");
  
  shift_float_to_complex(nz*ny*nx, (float *) datain);
  if(verbose>1) fprintf(stdout, "Rearranged data from float to complex\n");
  
  fftwf_execute(plan3d);
  if(verbose>1) fprintf(stdout, "Executed 3D FFT\n");
  
  /* Now scale the data since FFTW does not*/
  scale = 1.0/sqrt(1.0*nz*nx*ny);
  for(i=0;i<nz*ny*nx;i++)
    dataout[i] *= scale;
  
  /* Write that out to the file */
  result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)dataout);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Error in setting values.\n");
    return 0;
  }
  

  fftwf_destroy_plan(plan3d);   
  fftwf_free(datain); fftwf_free(dataout); 
  
  if(verbose)
    fprintf(stdout, "Completed full 3d FFT\n");
    
  return 1;
}

int full_3difft(mihandle_t filein, mihandle_t fileout, int flags)
{
  long mz,my,mx;
  long nz,ny,nx;
  long slice, x,y,z, i;
  fdr_complex *datain = NULL;
  fdr_complex *dataout = NULL;
  fftwf_plan plan3d;
  unsigned long start[3], count[3],nelems;
  float scale;
  
  int result;
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  
  // Dimensions are the same size, we can feel free to do the transforms!

  datain = (fdr_complex *)fftwf_malloc(nz*ny*nx*sizeof(fdr_complex));
  dataout = (fdr_complex *)fftwf_malloc(nz*ny*nx*sizeof(fdr_complex));
  
  plan3d = fftwf_plan_dft_3d(mz,my,mx, datain, dataout, FFTW_BACKWARD, FFTW_MEASURE);
  
  start[0] = (unsigned long)0; start[1] = start[2] = 0;
  count[0] = (unsigned long)mz;
  count[1] = (unsigned long)my;
  count[2] = (unsigned long)mx;

  result = miget_voxel_value_hyperslab(filein,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) datain);
  if(result == MI_ERROR) {
    fprintf(stderr, "Error getting data.\n");
    return 0;
  }
  else if(verbose>1) fprintf(stdout,"Retrieved temp3d data\n");
  
  fftwf_execute(plan3d);
  if(verbose>1) fprintf(stdout, "Executed 3D FFT\n");
  
  /* Now scale the data since FFTW does not*/
  scale = 1.0/sqrt(1.0*nz*nx*ny);
  for(i=0;i<nz*ny*nx;i++)
    dataout[i] *= scale;
  
  /* Write that out to the file */
  result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)dataout);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Error in setting values.\n");
    return 0;
  }

  fftwf_destroy_plan(plan3d);   
  fftwf_free(datain); fftwf_free(dataout); 
  
  if(verbose)
    fprintf(stdout, "Completed full 3d FFT\n");
    
  return 1;
}


/***************************************************************
 ***************************************************************
 ****************** SLICEWISE FUNCTIONS ************************
 ***************************************************************
 ***************************************************************/

int slicewise_print_data() {}

int slicewise_build_rolloff_filter(mihandle_t filter_file, double weight, double maxslope) 
{
  
  long n1, n2, n3, i,j,k,index1,index2;
  unsigned long start[3], count[3];
  float startk,starti,curk,curi,dk,di,slope;
  float *rolloff_filter;
  int result;
    
  result = get_minc_dimensions_from_handle(filter_file, &n3, &n2, &n1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "=>Got all the dimensions\n");
  
  /* First we build the 2D rolloff filter that will be replicated at each y-step */
  /* We then iterate over all y and copy the data in place */
  
  /* Allocate memory for the 2D rolloff filter */
  rolloff_filter = (float *)fftwf_malloc(n3*n1*sizeof(float)); 

  /* Build the rolloff filter */
  // Changed July 6 to account for even and odd length data sets
  dk = 1/(2*M_PI);
  di = 1/(2.0);
  if(n3%2)
    startk = ( - (n3-1.0) / 2.0) / (2.0*M_PI);
  else
    startk = ( - (  n3  ) / 2.0) / (2.0*M_PI);
  if(n1%2)
    starti = ( - (n1-1.0) / 2.0) / (2.0);
  else
    starti = ( - (  n1  ) / 2.0) / (2.0);  

  for(k=0,curk=startk;k<n3;k++,curk+=dk) {
    for(i=0,curi=starti;i<n1;i++,curi+=di) {
      index1 = k*n1 + i;
      slope = -curk/curi;
      /* Special case -- along the zero-axis */
      //if( (fabs(curk)<=1e-8) || fabs(curi)<=(1e-8)) {
      if(fabs(curk)<=(1e-8)) {
        rolloff_filter[index1] = 1.0;
        continue;
      }
      if(fabs(curi)<=(1e-8) && fabs(curk) < 0.5) {
        rolloff_filter[index1] = 1.0;
        continue;
      }
      /* Should be no data in these areas */
      if(fabs(slope) >= 1.0) {
        rolloff_filter[index1] = 0.0;
        continue;
      }
      /* Cosine rolloff from maxslope down to 1.0 */
      if(slope >= maxslope) {
        rolloff_filter[index1] = pow(cos(M_PI/2 * (fabs(slope)-maxslope) / (1.0-maxslope)), 2);
        continue;
      }
      /* Zero out the out of focus data beyond weight */
      if(slope <= -weight) {
        rolloff_filter[index1] = 0.0;
        continue;
      }
      /* Want to keep all data in the good quadrant*/
      if(slope >= 0) {
        rolloff_filter[index1] = 1.0;
        continue;
      }
      /* Everything else should be eliminated. We should be in rolloff region. */
      /* It is 1.0 at slope = 0.0, and 0.0 at slope = weight */
      rolloff_filter[index1] = pow(cos( M_PI/2.0 * slope/weight ),2);
    }

  }
  if(verbose>1) fprintf(stdout, "=>Constructed rolloff filter\n");
  /* Iterate over each y-step, setting the values then writing them out to file */
  for(j=0;j<n2;j++) {
    
    start[0] = 0; start[2] = 0; 
    start[1] = (unsigned long)j;
    count[0] = (unsigned long)n3;
    count[1] = 1;
    count[2] = (unsigned long)n1;
    
    result = miset_voxel_value_hyperslab(filter_file,MI_TYPE_FLOAT,start,count,(float *)rolloff_filter);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error in setting values.\n");
      return 0;
    }
    else if(verbose >1)
      fprintf(stdout, "Set values at y=%d.\n", j);
    
//    for(j=0;j<n2;j++) {
//      for(i=0;i<n1;i++) {
//        index1 = k*n2*n1 + j*n1 + i;
//        index2 = k*n1 + i;
//        outdata[index1] = rolloff[index2];
//      }
//    }
  } 
  
  result = miset_volume_valid_range(filter_file, 1.0,0.0);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Could not set volume range\n");
    return 0;
  } 

  fftwf_free(rolloff_filter);

  if(verbose>1) fprintf(stdout, "****** Returning from build_rolloff_filter ******\n");

  return 1;

}

int slicewise_build_wiener_filter(mihandle_t filein, mihandle_t fileout, long newelems) 
{
  long nelems,i,j,k, index;
  long n3, n2, n1, m3, m2, m1;
  int result;
  float noise;
  float *tmp;
  float *rad, *rcount;
  float radius, frac;
  float z,y,x;
  float scale;
  long iradius;
  unsigned long start[3], count[3];

  result = get_minc_dimensions_from_handle(filein, &n3, &n2, &n1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "=>Got dimensions from file 1\n");
  
  result = get_minc_dimensions_from_handle(fileout, &m3, &m2, &m1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "=>Got dimensions from file 2\n");
  
  if(m3 != n3 || m2 != n2 || m1 != n1)
  {
    fprintf(stdout, "=>Dimensions do not match up!");
    return 0;
  }

  rad = (float *)fftwf_malloc(newelems*sizeof(float));
  rcount = (float *)fftwf_malloc(newelems*sizeof(float));
  if(rad == NULL || count == NULL) {
    fprintf(stderr, "Could not allocate memory for radial power spectrum\n");
    return 0;
  }

  nelems = n3*n2*n1;

  tmp = (float *)fftw_malloc(3*3*3*sizeof(float));
  
  /* Let's get a mean value for the noise floor */
  /* For this we'll use a bit of data just off the center line */

  start[0] = n3/2-50; start[1] = start[2] = 0;
  count[0] = count[1] = count[2] = 3;
  
  // fprintf(stdout, "Selecting from %d, %d, %d along %d %d %d\n", start[0], start[1], start[2], count[0], count[1], count[2]);
  result = miget_voxel_value_hyperslab(filein, MI_TYPE_FLOAT, start, count, (float *)tmp);
  if(result == MI_ERROR) 
  { 
    fprintf(stdout, "Could not retrieve data\n"); fftwf_free(rad); fftwf_free(count); fftwf_free(tmp); return 0; 
  }
  if(verbose>1) fprintf(stdout, "Retrieved data to calculate noise\n");

  noise = 0;
  for(i=0;i<3*3*3;i++)
  {
    // fprintf(stdout, "at %d it's %2.2f\n", i, tmp[i]);
    noise += tmp[i];
  }
  noise /= (3.*3.*3.);
  if(verbose>1) fprintf(stdout, "Calculated noise floor of %f\n", noise);

  fftwf_free(tmp);
  
  /* Noise is now calculated, so let's begin doing the radially averaged power spectrum */
  scale = (newelems-2) / sqrt( (n3-1.)/2.*(n3-1)/2. + (n2-1.)/2.*(n2-1)/2. + (n1-1.)/2.*(n1-1.)/2. );
  tmp = fftwf_malloc(n2*n1*sizeof(float));
  /* Iterate over each z-slice in the power spectrum */
  for(k=0;k<n3;k++) {
    
    /* Retrieve the data from the file */
    start[0] = (unsigned long) k; start[1] = start[2] = 0;
    count[0] = 1;
    count[1] = (unsigned long)n2; count[2] = (unsigned long)n1;
    result = miget_voxel_value_hyperslab(filein, MI_TYPE_FLOAT, start, count, (float *)tmp);
    if(result == MI_ERROR) 
    { 
      fprintf(stdout, "Could not retrieve data from z=%d\n", k);
      fftwf_free(rad); fftwf_free(count); fftwf_free(tmp); return 0; 
    }
    if(verbose>1) fprintf(stdout, "Retrieved data from z=%d\n", k);
    
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
        index = j*n1 + i;
        radius = sqrt(z*z + y*y + x*x) * scale;
        iradius = (int)floor(radius);
        frac = radius - iradius;
        if(iradius < newelems-1) {
          rad[iradius] += (1-frac)*tmp[index];
          rcount[iradius] += (1-frac);
        }
        if(iradius+1 < newelems-1) {
          rad[iradius+1] += frac*tmp[index];
          rcount[iradius+1] += frac;
        }
      }
    }
  }
  
  fftwf_free(tmp);

  /* Correct for any zero values */
  for(i=0;i<newelems;i++) {
    if(rcount[i] == 0) 
      rad[i] = 0;
    else
      rad[i] = rad[i]/rcount[i];
  }

  /* Now we can build the filter itself */
  /* This is the radial filter, have to copy it to 3d later */
  for(i=0;i<newelems;i++) {
    /* This is the traditional Wiener filter */
    /* Our "signal" is just power spec of image minus noise */
    if(rad[i] != 0) {
      rad[i] = 1. /  (1. + noise / (rad[i]-noise) );
      if(rad[i] < 0)
      rad[i] = 0;
    }
    // fprintf(stdout, "The value at %d is %2.5f\n",  i, rad[i]);
  }
  // rad[newelems-1] = 1.0;

  /* So we have the 1-D radial Wiener filter.  Let's copy that into the 3D space */
  tmp = fftwf_malloc(n2*n1*sizeof(float));
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
        index = j*n1 + i;
        radius = sqrt(z*z + y*y + x*x) * scale;
        iradius = (int)floor(radius);
        frac = radius - iradius;
        tmp[index] = (1-frac)*rad[iradius] + frac*rad[iradius+1];
        if(tmp[index] < 0) {
          fprintf(stdout,"Bug?  frac=%f r=%f ir=%d rad=%f rad2=%f\n", frac, radius, iradius, rad[iradius],rad[iradius+1]);
          tmp[index] = 0;
        }
        if(tmp[index] > 1.0) {
          fprintf(stdout, "Bug?  frac=%f r=%f ir=%d rad=%f rad2=%f\n", frac, radius, iradius, rad[iradius],rad[iradius+1]);
          tmp[index] = 1;
        }
      }
    }
    start[0] = (unsigned long) k; start[1] = start[2] = 0;
    count[0] = 1;
    count[1] = (unsigned long)n2; count[2] = (unsigned long)n1;
    result = miset_voxel_value_hyperslab(fileout, MI_TYPE_FLOAT, start, count, (float *)tmp);
    if(result == MI_ERROR)
    {
      fftwf_free(tmp); fftwf_free(rcount); fftwf_free(rad);
      fprintf(stderr, "Could not set 3D Wiener filter data at z=%d\n", k);
      return 0;
    }
    if(verbose>1) fprintf(stdout, "Set values for 3D Wiener filter at z=%d\n", k);
  }
  
  result = miset_volume_valid_range(fileout, 1.0,0.0);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Could not set volume range\n");
    return 0;
  } 
  
  fftwf_free(rad); fftwf_free(rcount); fftwf_free(tmp);
  
  if(verbose>1) fprintf(stdout, "****** Calculated wiener filter ******\n");

  return 1;
}

int slicewise_build_bw_filter(mihandle_t filter_file, double bandlimit) 
{
  long n3, n2, n1, k,j,i, index;
  unsigned long start[3], count[3];
  float startk,dk, startj,dj, starti,di, curk, curj, curi;
  float bw_i, bw_j, scale;
  float *bw_filter;
  int result;

  // fprintf(stdout, "Received bandlimit of %2.2f\n", bandlimit2);
  
  result = get_minc_dimensions_from_handle(filter_file, &n3, &n2, &n1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "=>Got all the dimensions\n");

  /* Allocate memory for the 2D bw filter */
  bw_filter = (float *)fftwf_malloc(n2*n1*sizeof(float)); 

  // Changed July 6 to take odd and even numbered bins into account
  if(n3%2)
    startk = (-(n3-1.0)/2.0) / (2.0*M_PI);
  else
    startk = (-(  n3  )/2.0) / (2.0*M_PI);
  dk = 1/(2*M_PI);
  if(n2%2)
    startj = (-(n2-1.0)/2.0) / (2.0);
  else
    startj = (-(  n2  )/2.0) / (2.0);
  dj = 1/(2.0);
  if(n1%2)
    starti = (-(n1-1.0)/2.0) / (2.0);
  else
    starti = (-(  n1  )/2.0) / (2.0);
  di = 1/(2.0);

  // Find the right bandlimits
  bw_i = bandlimit*(n1-1)/4.0;
  bw_j = bandlimit*(n2-1)/4.0;
  
//  fprintf(stdout, "Starting at %2.2f, %2.2f, %2.2f with steps of %2.2f, %2.2f, %2.2f to end at %2.2f, %2.2f, %2.2f\n", 
//    startk, startj, starti, dk, di, dj, startk+n3*dk, startj+n2*dj, starti+n1*di);
//  fprintf(stdout, "Bandlimits %2.2f are along y: %2.4f and along x %2.4f", bandlimit2, bw_i, bw_j);
  for(k=0,curk=startk;k<n3;k++,curk+=dk) {
    for(j=0,curj=startj;j<n2;j++,curj+=dj) {
      for(i=0,curi=starti;i<n1;i++,curi+=di) {
        
        // We only need the index for the 2D file as we'll be writing it out at each z (or k)
        index = j*n1 + i;
    
        // At this point we want to check the bandlimit of the data
        // Three cases: 
        //     less than bandlimit*0.9 ==> scale=1
        //     less than bandlimit     ==> scale=cosine rolloff
        //     greater than bandlimit  ==> scale = 0

        // if(verbose>1) fprintf(stdout, "Bandlimit i is %f and j is %f at %f and %f\n", bw_i, bw_j, curi, curj);
        
        // fprintf(stdout, "Abs values are %2.2f, %2.2f, %2.2f\n", fabs(curi), fabs(curj), fabs(curk));
        if(fabs(curi) > bw_i || fabs(curj) > bw_j) {
          scale = 0.0;
        }
        else if(fabs(curi) < 0.9*bw_i && fabs(curj) < 0.9*bw_j) {
          // fprintf(stdout, "hello!\n");
          scale = 1.0;
        }
        
        else if( (fabs(curi) >= 0.9*bw_i && fabs(curi) <= bw_i) &&  (fabs(curj) >= 0.9*bw_j && fabs(curj) <= bw_j) ) {
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
        bw_filter[index] = scale;
        // fprintf(stdout, "Scale is %2.4f (%2.4f)\n", bw_filter[index], scale);
      }
    }
    
    start[0] = (unsigned long)k; start[1] = start[2] = 0;
    count[0] = 1; 
    count[1] = (unsigned long)n2;
    count[2] = (unsigned long)n1;
    
    // Write out the 2D data to file
    result = miset_voxel_value_hyperslab(filter_file,MI_TYPE_FLOAT,start,count,(float *)bw_filter);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error in setting values.\n");
      return 0;
    }
    if(verbose >1) fprintf(stdout, "Set values at z=%d.\n", k);
  }
  
  result = miset_volume_valid_range(filter_file, 1.0,0.0);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Could not set volume range\n");
    return 0;
  } 

  if(verbose) fprintf(stdout, "***** Completed build_bw_filter *****\n");

}

int slicewise_build_inverse_fdr_filter(mihandle_t psffile, mihandle_t filter_file, int flags) 
{
  long nelems;
  unsigned long psfindex1,psfindex2, start[3], count[3];
  long k,j,i, index;
  float startk,dk, startj,dj, starti,di, curk, curj, curi;
  float mindepth, maxdepth, ddepth, frac;
  long idepth;
  float depth, slope;
  double tmp1, tmp2, tmp3, denom1, denom2;
  float bw_i, bw_j, scale, tmp;
  long n3,n2,n1,psfn3,psfn2,psfn1;
  fdr_complex *psfdata, *filterdata;
  fdr_complex tmpv;
  int result;
  fdr_complex psf1, psf2;
  
  result = get_minc_dimensions_from_handle(psffile, &psfn3, &psfn2, &psfn1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "=>Retrieved dimensions from PSF file\n");
  
  result = get_minc_dimensions_from_handle(filter_file, &n3, &n2, &n1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "=>Retrieved dimensions from filter file\n");
  
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

  nelems = n3*n2*n1;

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
  
  ddepth = 2.0 / (psfn3-1);
  maxdepth = 1.0;
  mindepth = -1.0;
  fprintf(stdout, "Details are: startk %f dk %f startj %f dk %f starti %f di %f\n", startk, dk, startj,dj,starti,di);
  fprintf(stdout, "Depth details are: max %f min %f ddepth %f\n", maxdepth, mindepth, ddepth);
  // fprintf(stdout, "Bandlimit at: %f\n", bandlimit);
  fprintf(stdout, "PSF DIMS: Z=%d Y=%d X=%d\n", psfn3, psfn2, psfn1);
  fprintf(stdout, "SG  DIMS: Z=%d Y=%d X=%d\n", n3, n2, n1);

  /* Allocate memory */
  filterdata = (fdr_complex *)fftwf_malloc(n2*n1*sizeof(fdr_complex));
  psfdata = (fdr_complex *)fftwf_malloc(psfn3*psfn2*psfn1*sizeof(fdr_complex));
  start[0] = start[1] = start[2] = 0;
  count[0] = (unsigned long)psfn3; count[1] = (unsigned long)psfn2; count[2] = (unsigned long)psfn1;
  result = miget_voxel_value_hyperslab(psffile, MI_TYPE_FCOMPLEX, start, count, (fdr_complex *)psfdata);
  if(result == MI_ERROR)
  {
    fprintf(stdout, "Could not load the entire PSF file\n");
    fftwf_free(psfdata); fftwf_free(filterdata);
    return 0;
  }
  if(verbose>1) fprintf(stdout, "=>Loaded PSF data\n");
  
  /* Loop over all elements to get the appropriate information */
  /* Note this is iterating in frequency space, so k,j,i represent the coordinates denoted above */
  for(k=0,curk=startk;k<n3;k++,curk+=dk) {
    if(verbose) fprintf(stdout, "Now sampling Z=%d\n", k);

    for(j=0,curj=startj;j<n2;j++,curj+=dj) {
      for(i=0,curi=starti;i<n1;i++,curi+=di) {
        index = j*n1 + i;
        if(index>nelems) {
          fprintf(stdout, "WHOA!  index > nelems, this shouldn't be possible!\n");
          continue;
        }
        slope = -curk/curi;
        scale = 1.0;
        
        /* Just zero out everything that isn't "interesting */
        if(abs(slope) >= 1.0) {
          filterdata[index] = (fdr_complex)1.0;
          continue;
        }
        if(isinf(slope) || isnan(slope)) {
          filterdata[index] = (fdr_complex)1.0;
          continue;
        }
        
        depth = slope + 1.0;
        idepth = floor(depth/ddepth); 
        frac = depth/ddepth-idepth;
        count[0] = count[1] = count[2] = 1;
        psfindex1= idepth*psfn2*psfn1 + j*psfn1 + i;
        psfindex2 = (idepth+1)*psfn2*psfn1 + j*psfn1 + i;
        if(psfindex1 > psfn3*psfn2*psfn1)
        {
          fprintf(stderr, "Erm, wha?  psfindex1 > psfnelems?  %d vs %d (%f - %d depth) at %d,%d,%d\n", psfindex1, psfn3*psfn2*psfn1, depth, idepth, k,j,i);
        }
        if(psfindex2 > psfn3*psfn2*psfn1)
        {
          fprintf(stderr, "Erm, wha?  psfindex1 > psfnelems?  %d vs %d (%f - %d depth) at %d,%d,%d\n", psfindex2, psfn3*psfn2*psfn1, depth, idepth, k,j,i);
        }
        psf1 = (fdr_complex)psfdata[psfindex1];
        psf2 = (fdr_complex)psfdata[psfindex2];
        /*
        psfindex1[0] = (unsigned long) idepth;
        psfindex2[0] = (unsigned long) (psfindex1[0] + 1);
        psfindex1[1] =  psfindex2[1] = (unsigned long)j;
        psfindex1[2] = psfindex2[2] = (unsigned long)i;
        // fprintf(stdout, "Pulling from %d %d %d and %d %d %d\n", psfindex1[0], psfindex1[1], psfindex1[2], psfindex2[0], psfindex2[1], psfindex2[2]);
        result = miget_voxel_value_hyperslab(psffile, MI_TYPE_FCOMPLEX, psfindex1, count, &psf1);
        if(result == MI_ERROR)
        {
          fprintf(stdout, "Failed read of psf file.\n");
          return 0;
        }
        result = miget_voxel_value_hyperslab(psffile, MI_TYPE_FCOMPLEX, psfindex2, count, &psf2);
        if(result == MI_ERROR)
        {
          fprintf(stdout, "Failed read of psf file.\n");
          return 0;
        }
        */
        /*
          fprintf(stdout, "At coords Z=%d Y=%d X=%d we have curk=%f curj=%f, curi=%f\n", k, j, i, curk, curj, curi);
          
          fprintf(stdout, "Slope is %f at Z=%d, Y=%d, X=%d\n", slope, k, j, i);
          fprintf(stdout, "This corresponds to a depth of %f (index %d)\n", depth, idepth);
          fprintf(stdout, "And a psfindex of %d and %d\n", psfindex1, psfindex2);
        */
        /* First, let's interpolate to get the right frequency space value */
        tmpv = ((1-frac)*creal(psf1) + frac*creal(psf2)) + ((1-frac)*cimag(psf1) + frac*cimag(psf2))*I;
        tmp1 = creal(tmpv);
        tmp2 = cimag(tmpv);
        
        if(flags == INVERSE_FILTER) {
          /* 1 / (a+bi) = a/(a^2 + b^2) - bi/(a^2 + b^2) */
          if(!tmp1 && !tmp2) {
            filterdata[index] = (fdr_complex)0.0;
          }
          else {
            filterdata[index] = 1.0 / tmpv; // tmp1 / ( pow(tmp1,2)+pow(tmp2,2) ) * scale + (-tmp2 / ( pow(tmp1,2)+pow(tmp2,2) ) * scale)*I;
          }
          tmp3 = sqrt(pow(creal(filterdata[index]),2) + pow(cimag(filterdata[index]),2)); //sqrt(pow(cimag(filterdata[index]),2)+pow(creal(filterdata[index]),2));
          if(tmp3 < 0.0)
          {
            fprintf(stderr, "Why would this be less than zero?  %f (%f + %fi), %d,%d,%d\n", tmp3, creal(filterdata[index]), cimag(filterdata[index]), k,j,i);
          }
          if(scale == 1.0 &&  (float)tmp3 < 1.0) {
            fprintf(stderr, "BOOP BOOP ERROR!  fdrfilter mag == %f (%f +%fi) at value z,y,x=%d,%d,%d (%f,%f,%f).\n",tmp3, creal(filterdata[index]), cimag(filterdata[index]), k,j,i, curk, curj, curi);
            fprintf(stderr, "psfindex1=%d, psfindex2=%d, , depth = %f. idepth = %d, frac=%f.\n", 
                psfindex1, psfindex2, depth, idepth, frac);
            fprintf(stderr, "fdr[psfindex1]= %f + %f i, fdr[psfindex2] = %f + %f i\n", 
                creal(psf1), cimag(psf1), creal(psf2), cimag(psf2));
            fprintf(stderr, "Slope is %f, scale is %f\n", slope, scale);
            filterdata[index] = tmpv/tmp3; // creal(filterdata[index]) * 1./tmp1 +  cimag(filterdata[index]) * 1./tmp1 * I;
          
          }
          if(isinf(tmp1) || isinf(tmp2))
          {
            fprintf(stderr, "BOOP BOOP ERROR!  fdrfilter mag == %f (%f +%fi) at value z,y,x=%d,%d,%d (%f,%f,%f).\n",tmp3, creal(filterdata[index]), cimag(filterdata[index]), k,j,i, curk, curj, curi);
            fprintf(stderr, "psfindex1=%d, psfindex2=%d, , depth = %f. idepth = %d, frac=%f.\n", 
                psfindex1, psfindex2, depth, idepth, frac);
            fprintf(stderr, "fdr[psfindex1]= %f + %f i, fdr[psfindex2] = %f + %f i\n", 
                creal(psf1), cimag(psf1), creal(psf2), cimag(psf2));
            fprintf(stderr, "Slope is %f, scale is %f\n", slope, scale);           
            // fprintf(stderr, "Why would this be infinite?  %f (%f + %fi), %d,%d,%d\n", tmp3, creal(filterdata[index]), cimag(filterdata[index]), k,j,i);
            filterdata[index] = (fdr_complex)1.0;
          }
          if(isnan(tmp1) || isnan(tmp2))
          { 
            fprintf(stderr, "BOOP BOOP ERROR!  fdrfilter mag == %f (%f +%fi) at value z,y,x=%d,%d,%d (%f,%f,%f).\n",tmp3, creal(filterdata[index]), cimag(filterdata[index]), k,j,i, curk, curj, curi);
            fprintf(stderr, "psfindex1=%d, psfindex2=%d, , depth = %f. idepth = %d, frac=%f.\n", 
                psfindex1, psfindex2, depth, idepth, frac);
            fprintf(stderr, "fdr[psfindex1]= %f + %f i, fdr[psfindex2] = %f + %f i\n", 
                creal(psf1), cimag(psf1), creal(psf2), cimag(psf2));
            fprintf(stderr, "Slope is %f, scale is %f\n", slope, scale);
            // fprintf(stderr, "Why would this be nan?  %f (%f + %fi), %d,%d,%d\n", tmp3, creal(filterdata[index]), cimag(filterdata[index]), k,j,i);
            filterdata[index] = (fdr_complex)1.0;
          }
          /* Uncomment for testing purposes
          if(i==90 && abs(curj) == abs(dj*90))
          {
            fprintf(stdout, "Test: At %d %d %d we have %f + %f i (%f)\n", k,j,i,creal(filterdata[index]),cimag(filterdata[index]), (float)(tmp3));
          }
          */
        }
        else {
          filterdata[index] = scale * tmpv;
        }
      }
    }
    start[0] = (unsigned long)k;
    start[1] = start[2] = 0;
    count[0] = 1;
    count[1] = (unsigned long)n2; count[2] = (unsigned long)n1;
    
    result = miset_voxel_value_hyperslab(filter_file,MI_TYPE_FCOMPLEX,start,count, (mifcomplex_t *)filterdata);
    if(result == MI_ERROR) {
      fprintf(stderr, "Error writing data at z=%d.\n", k);
      return 0;
    }
    else if(verbose>1) fprintf(stdout,"=>Wrote filter data at z=%d\n", k);
    
  }
  
  //fprintf(stdout, "First value is %f + %f.i\n", fdpfilterdata[index].re, fdpfilterdata[index].im);
  if(verbose) fprintf(stdout,"Completed building the FDR filter.\n");
  
  /* Finally, we're successful */
  if(verbose) fprintf(stdout, "*********** Returning from build_inverse_fdr_filter **************\n");
  
  return 1;

}

int slicewise_limit_recovery_filter(mihandle_t filterfile, mihandle_t outfile, double fdrlimit) 
{
  long i, n3, n2, n1, m3, m2, m1,  z, nelems;
  float val1, val2; //, angle;
  int result;
  unsigned long count[3], start[3];
  fdr_complex *data;

  // fdrlimit = 5;
  fprintf(stdout, "Limiting with value %2.3f\n", fdrlimit);
  
  result = get_minc_dimensions_from_handle(filterfile, &n3, &n2, &n1);
  if(!result) { return 0; }
  if(verbose>1) fprintf(stdout, "=>Retrieved dimensions.\n");
  
  result = get_minc_dimensions_from_handle(filterfile, &m3, &m2, &m1);
  if(!result) { return 0; }
  if(verbose>1) fprintf(stdout, "=>Retrieved dimensions.\n");
  
  if(m3 != n3 || m2 != n2 || n1 != n1)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", m3,m2,m1, n3,n2,n1);
    return 0;
  }
 
  data = (fdr_complex *)fftw_malloc(n2*n1*sizeof(fdr_complex));
    
  nelems = n2*n1;
  for(z=0;z<n3;z++)
  {
    start[0] = (unsigned long)z;
    start[1] = start[2] = 0;
    count[0] = 1;
    count[1] = (unsigned long)n2; count[2] = (unsigned long)n1;
    result = miget_voxel_value_hyperslab(filterfile, MI_TYPE_FCOMPLEX, start, count, (fdr_complex *)data);
    if(result == MI_ERROR)
    {
      fprintf(stdout, "Could not retrieve data from the file\n");
      fftwf_free(data);
      return 0;
    }
    
    for(i=0;i<nelems;i++) {
      val1 = fabs(data[i]);
  
      if(val1 > fdrlimit/2) {
        val2 = fdrlimit - fdrlimit/2*exp(-(val1-fdrlimit/2)/fdrlimit);
        // fprintf(stdout, "val1 is: %2.2f val2 is %3.3f limit is %2.2f\n", val1, val2, fdrlimit);
        //fprintf(stdout, "Converting %2.2f + %2.2fi", creal(data[i]), cimag(data[i]));
        data[i] *= val2/val1;
        // fprintf(stdout, "to %2.2f + %2.2fi\n", creal(data[i]), cimag(data[i]));
      }
    }
    result = miset_voxel_value_hyperslab(outfile, MI_TYPE_FCOMPLEX, start, count, (mifcomplex_t *)data);
    if(result == MI_ERROR)
    {
      fprintf(stdout, "Could not set data in the file\n");
      fftwf_free(data);
      return 0;
    } 
    if(verbose>1) fprintf(stdout, "Set the data at z=%d\n", z); 
  }
  
  if(verbose>1) fprintf(stdout, "******** Completed limit_recovery_filter *********\n");
  return 1;
}

int slicewise_psfstack_fft(mihandle_t filein, mihandle_t fileout) 
{
  long mz,my,mx;
  long nz,ny,nx;
  long index1, index2, j, i;
  long slice, x,y,z;
  fdr_complex *temp3d = NULL;
  fdr_complex *datain2d = NULL;
  fdr_complex *dataout2d = NULL;
  fdr_complex *datain1d = NULL;
  fdr_complex *dataout1d = NULL;
  fdr_complex tmp;
  fftwf_plan plan2d, plan1d;
  unsigned long start[3], count[3],nelems;
  float scale;
  
  int result;
  
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
  
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
  
  if(mz != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }
  
  // Dimensions are the same size, we can feel free to do the transforms!

  temp3d    = (fdr_complex *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(fdr_complex)); 
  datain2d  = (fdr_complex *)fftwf_malloc(mx*my*sizeof(fdr_complex));
  dataout2d = (fdr_complex *)fftwf_malloc(nx*ny*sizeof(fdr_complex));
 
  plan2d = fftwf_plan_dft_2d(my,mx, (fdr_complex *)datain2d, (fdr_complex *)dataout2d, FFTW_FORWARD, FFTW_MEASURE);
  
  /* 
   * First we do a series of 2D FFTs along the z-direction, and store the
   * values in fileout
   */
  for(z=0;z<nz;z++)
  {
    /* 
     * Retrieve the data from the file 
     */
    start[0] = (unsigned long)z; start[1] = start[2] = 0;
    count[0] = 1;
    count[1] = (unsigned long)my;
    count[2] = (unsigned long)mx;
    /*
     * Assume that the filein is real, load it in, then shift to complex
     */
    result = miget_voxel_value_hyperslab(filein,MI_TYPE_FLOAT,start,count, (float *) datain2d);
    if(result == MI_ERROR) {
      fprintf(stderr, "Error getting data.\n");
      return 0;
    }
    else if(verbose>1) fprintf(stdout,"Retrieved temp3d data\n");
    
    /* Convert to complex */
    shift_float_to_complex(ny*nx, (float *) datain2d);
    if(verbose>1) fprintf(stdout, "Rearranged data from float to complex\n");
    
    scale = 1/sqrt(1.0*ny*nx);
//    for(y=0;y<ny;y++)
//    {
//      for(x=0;x<nx;x++)
//      {
//        datain2d[y*nx+x] = (float complex)temp3d[ny*nx+y*nx+x];
//      }
//    }
  
    /*
    for(y=0;y<ny*nx;y++)
    {
      fprintf(stdout, "Shifted result is (%d/%d): %f %f\n", y, ny*nx, creal(datain2d[y]), cimag(datain2d[y]));
    }
    */
    /* Perform the 2D FFT */
    fftwf_execute(plan2d);
    if(verbose>1) fprintf(stdout, "Executed 2D FFT at z=%d\n", z);
    /*
    for(y=0;y<ny*nx;y++)
    {
      fprintf(stdout, "Shifted result is (%d/%d): %f %f\n", y, ny*nx, creal(dataout2d[y]), cimag(dataout2d[y]));
    }
    */
    
    // 2-D FFT shift
    for(y=0;y<ny/2;y++)
    {
      for(x=0;x<nx/2;x++)
      {
        // swap quadrant 1 and 3
        index1 = y*nx + x;
        index2 = (y+ny/2)*nx + x+nx/2;
        tmp = dataout2d[index1];
        dataout2d[index1] = dataout2d[index2] * scale;
        dataout2d[index2] = tmp * scale;
        // swap quadrant 2 and 4
        index1 = y*nx + (x+nx/2);
        index2 = (y+ny/2)*nx + x;
        tmp = dataout2d[index1];
        dataout2d[index1] = dataout2d[index2] * scale;
        dataout2d[index2] = tmp * scale;
      }
    }
    // if(verbose>1) fprintf(stdout, "Shifted data.\n");
    
    /* Write that out to the file */
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)dataout2d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error in setting values.\n");
      return 0;
    }
    else if(verbose >1)
      fprintf(stdout, "Set values at z=%d.\n", z);
    
  }
  
  fftwf_destroy_plan(plan2d);   
  fftwf_free(temp3d); 
  fftwf_free(datain2d); fftwf_free(dataout2d); 
  
  if(verbose)
    fprintf(stdout, "Completed PSF stack of 2D FFTs\n");

}

int slicewise_normalize_complexpsfstack() {}

int slicewise_normalize_psfstack(mihandle_t filein, mihandle_t fileout) 
{
  // Changed Mar 17 to normalize each slice to 1
  long i, j, k;
  unsigned long start[3], count[3];
  long n3, n2, n1, m3, m2, m1;
  long index, index1, index2;
  float tempsum;
  int result;
  float *data;
  float dmax, dmin;

  dmax = -1e100;
  dmin = 1e100;
  result = get_minc_dimensions_from_handle(filein, &n3, &n2, &n1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "=>Retrieved minc dimensions from handle 1.\n");

  result = get_minc_dimensions_from_handle(fileout, &m3, &m2, &m1);
  if(!result) return 0;
  if(verbose>1) fprintf(stdout, "=>Retrieved minc dimensions from handle 2.\n");
  
  if(m3 != n3 || m2 != m2 || m1 != n1)
  {
    fprintf(stdout, "Dimensions do not match up!  exiting ... \n");
    return 0;
  }
  
  data = (float *) fftwf_malloc(n2*n1*sizeof(float));
  
  /* Iterate over depth point */
  for(k=0;k<n3;k++) {
    start[0] = (unsigned long)k;
    start[1] = start[2] = 0;
    count[0] = 1;
    count[1] = (unsigned long) n2; count[2] = (unsigned long)n1;
    
    result = miget_voxel_value_hyperslab(filein,MI_TYPE_FLOAT,start,count, (float *) data);
    if(result == MI_ERROR)
    {
      fftwf_free(data);
      fprintf(stdout, "Could not read data at z=%d\n", k);
      return 0;
    }
    if(verbose>1) fprintf(stdout, "Read values at z=%d\n", k);
    tempsum = 0;
    for(j=0;j<n2;j++) {
      for(i=0;i<n1;i++) {
        index = j*n1 + i;
        tempsum += data[index];
      }
    }
    for(j=0;j<n2;j++) {
      for(i=0;i<n1;i++) {
        index = j*n1 + i;
        data[index] /= tempsum;
        if(data[index] < dmin)
          dmin = data[index];
        if(data[index] > dmax)
          dmax = data[index];
      }
    }
   
    result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FLOAT,start,count, (float *) data);
    if(result == MI_ERROR)
    {
      fftwf_free(data);
      fprintf(stdout, "Could not read data at z=%d\n", k);
      return 0;
    }
    if(verbose>1) fprintf(stdout, "Wrote values at z=%d\n", k);
  }

  result = miset_volume_valid_range(fileout, dmax, dmin);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Could not set volume valid range\n");
    return 0;
  }
  result = miset_volume_range(fileout, dmax, dmin);
  if(result == MI_ERROR)
  {
    fprintf(stderr, "Could not set volume range\n");
    return 0;
  } 

  if(verbose>1) fprintf(stdout,"Completed PSF stack normalization\n");

  /* Finally, we're successful */
  return 1;


}

int slicewise_3dfft(mihandle_t filein, mihandle_t fileout, int flags)
{
	long mz,my,mx;
	long nz,ny,nx;
	long slice, x,y,z;
  fdr_complex *temp3d = NULL;
	fdr_complex *datain2d = NULL;
	fdr_complex *dataout2d = NULL;
	fdr_complex *datain1d = NULL;
	fdr_complex *dataout1d = NULL;
	fftwf_plan plan2d, plan1d;
	unsigned long start[3], count[3],nelems;
	float scale;
	
	int result;
  
	result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
	if(!result) return 0;
	
	result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
	if(!result) return 0;
	
	if(mz != nz || my != ny || nx != nx)
	{
		fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
		return 0;
	}
	
	// Dimensions are the same size, we can feel free to do the transforms!

  temp3d    = (fdr_complex *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(fdr_complex));	
	datain2d  = (fdr_complex *)fftwf_malloc(mx*my*sizeof(fdr_complex));
  datain1d  = (fdr_complex *)fftwf_malloc(mz*sizeof(fdr_complex));
	dataout2d = (fdr_complex *)fftwf_malloc(nx*ny*sizeof(fdr_complex));
  dataout1d = (fdr_complex *)fftwf_malloc(nz*sizeof(fdr_complex));
 
	plan2d = fftwf_plan_dft_2d(my,mx, datain2d, dataout2d, FFTW_FORWARD, FFTW_MEASURE);
  plan1d = fftwf_plan_dft_1d(mz, datain1d, dataout1d, FFTW_FORWARD, FFTW_PATIENT);
  
	
	/* 
	 * First we do a series of 2D FFTs along the z-direction, and store the
	 * values in fileout
	 */
	for(z=0;z<nz;z+=NUMSLICES)
	{
		/* 
		 * Retrieve the data from the file 
		 */
	  start[0] = (unsigned long)z; start[1] = start[2] = 0;
  	count[0] = NUMSLICES;
	  count[1] = (unsigned long)my;
  	count[2] = (unsigned long)mx;
  	/*
  	 * Assume that the filein is real, load it in, then shift to complex
  	 */
		result = miget_voxel_value_hyperslab(filein,MI_TYPE_FLOAT,start,count, (float *) temp3d);
	  if(result == MI_ERROR) {
			fprintf(stderr, "Error getting data.\n");
			return 0;
  	}
  	else if(verbose>1) fprintf(stdout,"Retrieved temp3d data\n");
    
    /* Convert to complex */
    shift_float_to_complex(NUMSLICES*ny*nx, (float *) temp3d);
    if(verbose>1) fprintf(stdout, "Rearranged data from float to complex\n");
     
    for(slice=0;slice<NUMSLICES;slice++)
    {
      for(y=0;y<ny;y++)
      {
        for(x=0;x<nx;x++)
        {
          datain2d[y*nx+x] = (float complex)temp3d[slice*ny*nx+y*nx+x];
        }
      }
  
      /*
  		for(y=0;y<ny*nx;y++)
  		{
  			fprintf(stdout, "Shifted result is (%d/%d): %f %f\n", y, ny*nx, creal(datain2d[y]), cimag(datain2d[y]));
  		}
  		*/
  		/* Perform the 2D FFT */
  		fftwf_execute(plan2d);
  		if(verbose>1) fprintf(stdout, "Executed 2D FFT at z=%d\n", z+slice);
  		/*
  		for(y=0;y<ny*nx;y++)
  		{
  			fprintf(stdout, "Shifted result is (%d/%d): %f %f\n", y, ny*nx, creal(dataout2d[y]), cimag(dataout2d[y]));
  		}
      */
      for(y=0;y<ny;y++)
      {
        for(x=0;x<nx;x++)
        {
          temp3d[slice*ny*nx+y*nx+x] = (float complex)dataout2d[y*nx+x];
        }
      }
    }
		/* Write that out to the file */
		result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)temp3d);
		if(result == MI_ERROR)
		{
			fprintf(stderr, "Error in setting values.\n");
			return 0;
		}
		else if(verbose >1)
			fprintf(stdout, "Set values in the block beginning at z=%d.\n", z);
		
	}
	
  fftwf_destroy_plan(plan2d);   
  fftwf_free(temp3d); 
  fftwf_free(datain2d); fftwf_free(dataout2d); 
  
  /*
   * Reinitialize the temp3d buffer with different # elements
   */
  temp3d    = (fdr_complex *)fftwf_malloc(NUMSLICES*nz*nx*sizeof(fdr_complex));
	scale = 1.0/sqrt(1.0*nz*nx*ny);
  
	for(y=0;y<ny;y+=NUMSLICES)
	{
		/*
		 * Load up the 2D slices here in order to avoid loading up nx number of
		 * 1D slices below
		 */
		start[0] = 0;
		start[1] = (unsigned long)y;
		start[2] = 0;
		count[0] = (unsigned long)nz;
		count[1] = NUMSLICES;
    count[2] = (unsigned long)nx;
		result = miget_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(fdr_complex *) temp3d);
		if(result == MI_ERROR)
		{
			fprintf(stderr, "Error opening 1D values.\n");
			return 0;
		}
    if(verbose>1) fprintf(stdout, "Read in 2D slices\n");
		
    for(slice=0;slice<NUMSLICES;slice++)
    {
  		if(verbose>1) fprintf(stdout,"Now beginning 1D FFTs along y=%d\n", y+slice);
  		for(x=0;x<nx;x++)
  		{
  			/*
  			 * Now initialize the 1d array
  			 */
        for(z=0;z<nz;z++)
        {
          datain1d[z] = temp3d[z*NUMSLICES*nx+slice*nx+x];
        }
        // if(verbose>1) fprintf(stdout, "Initialized the 1D array\n");
        
  			fftwf_execute(plan1d);
        // if(verbose>1) fprintf(stdout, "Executed FFT\n");
  			/*
  			 * Need to perform intensity scaling
  			 */
        for(z=0;z<nz;z++)
        {
          temp3d[z*NUMSLICES*nx+slice*nx+x] = dataout1d[z] * scale;
        }
        // if(verbose>1) fprintf(stdout, "Scaled data\n");
      }
  			
  	}
    if(verbose>1) fprintf(stdout, "About to write 2D collection of 1D values\n");
    result = miset_voxel_value_hyperslab(fileout, MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *) temp3d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error writing 2D collection of 1D values.\n");
      return 0;
    }
		if(verbose>1) fprintf(stdout, "Completed 1D FFTs along y=%d\n", y);
	}
	
  fftwf_destroy_plan(plan1d);
	fftwf_free(temp3d); fftwf_free(datain1d); fftwf_free(dataout1d);
	
	if(verbose)
		fprintf(stdout, "Completed slicewise 3d FFT\n");
		
	return 1;
	
}

int slicewise_3difft(mihandle_t filein, mihandle_t fileout, int flags)
{
	long mz,my,mx;
	long nz,ny,nx;
	long x,y,z, slice;
  
  fdr_complex *temp3d = NULL;
	fdr_complex *datain2d = NULL;
	fdr_complex *dataout2d = NULL;
	fdr_complex *datain1d = NULL;
	fdr_complex *dataout1d = NULL;
	fftwf_plan plan2d, plan1d;
	unsigned long start[3], count[3],nelems;
	float scale;
	
	int result;
	
	result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
	if(!result) return 0;
	
	result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
	if(!result) return 0;
	
	if(verbose>1) fprintf(stdout, "Got all the dimensions\n");
	
	if(mz != nz || my != ny || nx != nx)
	{
		fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
		return 0;
	}
	
	// Dimensions are the same size, we can feel free to do the transforms!

  temp3d    = (fdr_complex *)fftwf_malloc(NUMSLICES*ny*nx*sizeof(fdr_complex));	
	datain2d  = (fdr_complex *)fftwf_malloc(mx*my*sizeof(fdr_complex));
	dataout2d = (fdr_complex *)fftwf_malloc(nx*ny*sizeof(fdr_complex));
	datain1d  = (fdr_complex *)fftwf_malloc(mz*sizeof(fdr_complex));
	dataout1d = (fdr_complex *)fftwf_malloc(nz*sizeof(fdr_complex));
	
	plan2d = fftwf_plan_dft_2d(my,mx, datain2d, dataout2d, FFTW_BACKWARD, FFTW_MEASURE);
	plan1d = fftwf_plan_dft_1d(mz, datain1d, dataout1d, FFTW_BACKWARD, FFTW_MEASURE);
	
	if(verbose>1) fprintf(stdout, "Created plans.\n");
	
	/* 
	 * First we do a series of 2D FFTs along the z-direction, and store the
	 * values in fileout
	 */
	for(z=0;z<nz;z+=NUMSLICES)
	{
		/* 
		 * Retrieve the data from the file 
		 */
	  start[0] = (unsigned long)z; start[1] = start[2] = 0;
  	count[0] = NUMSLICES;
	  count[1] = (unsigned long)my;
  	count[2] = (unsigned long)mx;
  	/*
  	 * Assume that the filein is real, load it in, then shift to complex
  	 */
		result = miget_voxel_value_hyperslab(filein,MI_TYPE_FCOMPLEX,start,count, (fdr_complex *) temp3d);
	  if(result == MI_ERROR) {
			fprintf(stderr, "Error getting data.\n");
			return 0;
  	}
  	else if(verbose>1) fprintf(stdout,"Retrieved data\n");

    for(slice=0;slice<NUMSLICES;slice++)
    {
      for(y=0;y<ny;y++)
      {      
        for(x=0;x<nx;x++)
        {
          datain2d[y*nx+x] = (fdr_complex) temp3d[slice*ny*nx+y*nx+x];
        }
      } // end setting 2d array in
      
    	/* Perform the 2D FFT */
    	fftwf_execute(plan2d);
    	if(verbose>1) fprintf(stdout, "Executed 2D IFFT at z=%d\n", z+slice);
      
      for(y=0;y<ny;y++)
      {      
        for(x=0;x<nx;x++)
        {
          temp3d[slice*ny*nx+y*nx+x] = (fdr_complex) dataout2d[y*nx+x];
        }
      } // end setting 2d array out
    } // end iterating over slices
      
		/* Write that out to the file */
		result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *)temp3d);
		if(result == MI_ERROR)
		{
			fprintf(stderr, "Error in setting values.\n");
			return 0;
		}
		else if(verbose >1)
			fprintf(stdout, "Set values for the block beginning with z=%d.\n", z);
		
	} // end iterating over z
  
  fftwf_destroy_plan(plan2d);
  fftwf_free(datain2d); fftwf_free(dataout2d); fftwf_free(temp3d);
  
  /* 
   * Reinitialize temp3d with different dimensions
   */
  temp3d    = (fdr_complex *)fftwf_malloc(NUMSLICES*nz*nx*sizeof(fdr_complex));  
	scale = 1.0/sqrt(1.0*nz*ny*nx);
  
	
	for(y=0;y<ny;y+=NUMSLICES)
	{
    /*
	   * Load up the 2D slices here in order to avoid loading up nx number of
	   * 1D slices below
	   */
	  start[0] = 0;
	  start[1] = (unsigned long)y;
	  start[2] = 0;
	  count[0] = (unsigned long)nz;
	  count[1] = NUMSLICES;
	  count[2] = (unsigned long)nx;
	  result = miget_voxel_value_hyperslab(fileout,MI_TYPE_FCOMPLEX,start,count,(fdr_complex *) temp3d);
	  if(result == MI_ERROR)
	  {
	    fprintf(stderr, "Error opening 1D values.\n");
	    return 0;
	  }
    
	  for(slice=0;slice<NUMSLICES;slice++)
    {
      if(verbose>1) fprintf(stdout,"Now beginning 1D FFTs along y=%d\n", y+slice);        
  		for(x=0;x<nx;x++)
  		{
  			/*
  			 * Set up the data for the transform from the already-read in 2D slice
  			 */
        for(z=0;z<nz;z++)
          datain1d[z] = temp3d[z*NUMSLICES*nx+slice*nx+x];
  			
  			fftwf_execute(plan1d);
  
  			/*
  			 * Need to perform intensity scaling and set back into 3d buffer
  			 */
        for(z=0;z<nz;z++)
        {
           temp3d[z*NUMSLICES*nx+slice*nx+x] = dataout1d[z] * scale;
        }
  		}// end looping over x
    } // end looping over slice
    
    result = miset_voxel_value_hyperslab(fileout, MI_TYPE_FCOMPLEX,start,count,(mifcomplex_t *) temp3d);
    if(result == MI_ERROR)
    {
      fprintf(stderr, "Error writing 3D collection of 1D values.\n");
      return 0;
    }
    if(verbose>1) fprintf(stdout, "Wrote 3D collection of 1D values beginning at y=%d.\n", y);
	}
	
  fftwf_destroy_plan(plan1d);
	fftwf_free(datain1d); fftwf_free(dataout1d); fftwf_free(temp3d);
	
	if(verbose)
		fprintf(stdout, "Completed slicewise 3d IFFT\n");
		
	return 1;
	
}

int slicewise_antialias_data(mihandle_t filein, mihandle_t fileout, long cutoff) 
{
  
  long mz,my,mx;
  long nz,ny,nx;
  long slice, x,y,z;
  fdr_complex *temp3din = NULL;
  fdr_complex *temp3dout = NULL;
  unsigned long startin[3], startout[3], count[3];
  
  int result;
    
  result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
  if(!result) return 0;
    
  result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
  if(!result) return 0;
    
  if(mz*3 != nz || my != ny || nx != nx)
  {
    fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
    return 0;
  }

  temp3din  = (fdr_complex *)fftwf_malloc(NUMSLICES*my*mx*sizeof(fdr_complex));
  temp3dout = (fdr_complex *)fftwf_malloc(NUMSLICES*my*mx*sizeof(fdr_complex));

  for(z=0;z<nz;z+=NUMSLICES)
  {
    if(z<cutoff || z > nz-cutoff)
    {
      fprintf(stdout, "Setting z= %d-%d to be zero.\n", z, z+NUMSLICES);
      /*
      * Don't need to read the file in -- these will be zeros in the output
      */
      startout[0] = (unsigned long)z; startout[1] = startout[2] = 0;
      count[0] = NUMSLICES;
      count[1] = (unsigned long)my;
      count[2] = (unsigned long)mx;  
      for(slice=0;slice<NUMSLICES;slice++)
      {
        for(y=0;y<ny;y++)
        {
          for(x=0;x<nx;x++)
          {
            temp3dout[slice*ny*nx+y*nx+x] = 0;
          }
        }
      } // end slice loop
      result = miset_voxel_value_hyperslab(fileout, MI_TYPE_FCOMPLEX,startout,count, (mifcomplex_t*) temp3dout);
      if(result == MI_ERROR) {
        fprintf(stderr, "Error setting data.\n");
        return 0;
      }
    } // end z loop
    else
    {
      fprintf(stdout, "Copying z=%d-%d from the input.\n", z, z+NUMSLICES);
      if(z < mz)
      {
        startin[0] = (unsigned long)z; startin[1] = startin[2] = 0;
        startout[0] = (unsigned long)z; startin[1] = startin[2] = 0;
        count[0] = NUMSLICES;
        count[1] = (unsigned long)my;
        count[2] = (unsigned long)mx;
      }
      else if(z < mz*2)
      {
        startin[0] = (unsigned long)z-mz; startin[1] = startin[2] = 0;
        startout[0] = (unsigned long)z; startin[1] = startin[2] = 0;
        count[0] = NUMSLICES;
        count[1] = (unsigned long)my;
        count[2] = (unsigned long)mx;
      }
      else
      {
        startin[0] = (unsigned long)z-mz*2; startin[1] = startin[2] = 0;
        startout[0] = (unsigned long)z; startin[1] = startin[2] = 0;
        count[1] = (unsigned long)my;
        count[2] = (unsigned long)mx;
      }
      /*
      * Assume that the filein is complex so read it in
      */    
      fprintf(stdout, "Collecting data from %d %d %d -- %d %d %d\n", startin[0], startin[1], startin[2],
              count[0], count[1], count[2]);
      result = miget_voxel_value_hyperslab(filein,MI_TYPE_FCOMPLEX,startin,count, (fdr_complex*) temp3din);
      if(result == MI_ERROR) {
        fprintf(stderr, "Error getting data.\n");
        return 0;
      }
      for(slice=0;slice<NUMSLICES;slice++)
      {   
        for(y=0;y<ny;y++)
        {
          for(x=0;x<nx;x++)
          {
            temp3dout[slice*ny*nx+y*nx+x] = temp3din[slice*ny*nx+y*nx+x];
          }
        }
      } // end slice loop
      result = miset_voxel_value_hyperslab(fileout, MI_TYPE_FCOMPLEX,startout,count, (mifcomplex_t*) temp3dout);
      if(result == MI_ERROR) {
        fprintf(stderr, "Error setting data.\n");
        return 0;
      }  
    } // end z loop
  }
      
    
  fftwf_free(temp3din);
  fftwf_free(temp3dout);
  return 1;
  
}

int slicewise_calculate_magnitude(mihandle_t filein, mihandle_t fileout, int flags) {

	long nz,ny,nx;
	long mz,my,mx;
	unsigned long start[3], count[3];
	long z,i;
	fdr_complex *tmp;
	float *tmp2;
	float x;
	fdr_complex *data;
	double dmin,dmax;

	int result;
	
	dmin = 1e10; dmax = 0;
	result = get_minc_dimensions_from_handle(filein, &mz, &my, &mx);
	if(!result) return 0;
	
	result = get_minc_dimensions_from_handle(fileout, &nz, &ny, &nx);
	if(!result) return 0;
	
	if(verbose>1) fprintf(stdout, "Got all the dimensions\n");
	
	if(mz != nz || my != ny || nx != nx)
	{
		fprintf(stderr, "Dimensions do not match (in: [%d,%d,%d] out: [%d,%d,%d])\n", mz,my,mx, nz,ny,nx);
		return 0;
	}
	
	data = (fdr_complex *)fftwf_malloc(nx*ny*sizeof(fdr_complex));
	
	for(z=0;z<nz;z++)
	{
	/* 
	 * First we do a series of 2D FFTs along the z-direction, and store the
	 * values in fileout
	 */
	 /* 
		 * Retrieve the data from the file 
		 */
	  start[0] = (unsigned long)z; start[1] = start[2] = 0;
  	count[0] = 1;
	  count[1] = (unsigned long)my;
  	count[2] = (unsigned long)mx;
		result = miget_voxel_value_hyperslab(filein,MI_TYPE_FCOMPLEX,start,count,(fdr_complex *) data);
	  if(result == MI_ERROR) {
			fprintf(stderr, "Error getting data.\n");
			return 0;
  	}
   /*
  	else if(verbose>1) 
  		fprintf(stdout,"Retrieved data from z=%d\n", z);
   */
  	/* 
  	 * We're going to do the magnitude in-place
  	 */	
  	tmp = (fdr_complex *)data;
  	tmp2 = (float *)data;
	  for(i=0;i<nx*ny;i++) {
      /* calculate power spectrum */
      if(flags)
        x = fabs(creal(tmp[i])*creal(tmp[i]) + cimag(tmp[i])*cimag(tmp[i]));
      /* calculate magnitude */
      else
  			x = (float)fabs(tmp[i]);
			tmp2[i] = x;
      /*
      if(verbose>1)
        fprintf(stdout, "The value at index=%d is %f from the value %f + %f i \n", i, tmp2[i], creal(tmp[i]), cimag(tmp[i]));
      */
			if(tmp2[i] < dmin)
				dmin = tmp2[i];
			if(tmp2[i] > dmax)
				dmax = tmp2[i];
			if(isinf(tmp2[i])) 
	  		fprintf(stdout, "Infinity calculated at index %d, from %f %f\n", i, creal(tmp[i]), creal(tmp[i]));
  	}
  		
		/* Write that out to the file */
		result = miset_voxel_value_hyperslab(fileout,MI_TYPE_FLOAT,start,count,(float *)data);
		if(result == MI_ERROR)
		{
			fprintf(stderr, "Error in setting values.\n");
			return 0;
		}
		else if(verbose >1)
			fprintf(stdout, "Set values at z=%d.\n", z);
		
	}
		
	if(verbose>1) fprintf(stdout, "Min: %f max: %f\n", dmin, dmax);
	fftwf_free(data);
	
	/* Set the range on the file */

	result = miset_volume_valid_range(fileout, dmax, dmin);
	if(result == MI_ERROR)
	{
		fprintf(stderr, "Could not set volume valid range\n");
		return 0;
	}
	result = miset_volume_range(fileout, dmax, dmin);
	if(result == MI_ERROR)
	{
		fprintf(stderr, "Could not set volume range\n");
		return 0;
	}	
	
	if(verbose)
		fprintf(stdout, "Completed slicewise magnitude\n");
		
	return 1;

}

int main(int argc, char *argv[]) {
	fprintf(stdout, "You shouldn't run this directly.\n");
	return 1;
}
