#include "fdr.h"

int isverbose=1;
int stepwiseoutput=1;

int main() {

  int result;
  long n3, n2, n1, nelems;
  long psfn3, psfn2, psfn1;
  fdr_complex *sgdata = NULL;
  fdr_complex *sgfftdata = NULL;
  fdr_complex *psfdata = NULL;
  fdr_complex *psffftdata = NULL;
  fdr_complex *fdpfilterdata = NULL;
  float        *noisedampendata = NULL;

  char *sgfile = "/projects/mice/jwalls/data/opt/fdp_simulation/3dtest/exp8/dataminusdark.mnc";
  char *psffile = "/projects/mice/jwalls/data/opt/fdp_simulation/3dtest/exp8/psf/PSF.mnc";
  char *finalfile = "/projects/mice/jwalls/data/opt/fdp_simulation/3dtest/exp8/finalback.mnc";

  if(isverbose) fprintf(stdout, "Beginning program ...\n");

  /* Get the sinogram data */
  result = open_minc_file(sgfile, &n3, &n2, &n1, &sgdata, REAL_AS_COMPLEX);
  if(!result) return 0;
  if(isverbose) fprintf(stdout, "--> Opened sinogram data\n");

  /* Get data of PSF stack */
  result = open_minc_file(psffile, &psfn3, &psfn2, &psfn1, &psfdata, REAL_AS_COMPLEX);
  if(!result) { free(sgfftdata); return 0; }
  if(isverbose) fprintf(stdout, "--> Opened psf stack\n");

  /* Let's check the dimensions */
  if(psfn2 != n2 || psfn1 != n1) {
	fprintf(stderr, "The X,Y dimensions do not match up.\n");
	fprintf(stderr,"SG: Y=%d, X=%d, PSF: Y=%d, X=%d\n", n2, n1, psfn2, psfn1);
	free(psfdata); free(sgdata);
  }
  if(isverbose) fprintf(stdout, "--> Dimensions match up, that's good.\n");
  
  nelems = n3*n2*n1;

  /* Allocate all the memory we need up front */
  fdpfilterdata = (fdr_complex *)calloc(nelems,sizeof(fdr_complex));
  noisedampendata = (float *)calloc(nelems, sizeof(float));
  if(fdpfilterdata==NULL || noisedampendata==NULL) {
	fprintf(stderr, "Could not allocate memory!\n");
	free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata);
	return 0;
  }
  if(isverbose) fprintf(stdout, "--> Allocated memory\n");

  /* Do the forward 3DFFT */
  result = sg_ccfft3d(n3,n2,n1,sgdata,FORWARD_FFT); /* Note this is an in-place transform! */
  if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(isverbose) fprintf(stdout, "--> Completed 3DFFT\n");
  sgfftdata = sgdata; /* Update pointer for naming convention's sake */
  if(stepwiseoutput) {
	result = write_minc_file("sgfft.mnc", n3, n2, n1, sgfftdata, COMPLEX_AS_COMPLEX);
	if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  }
  if(isverbose) fprintf(stdout, "--> Wrote file\n");

  /* Do the FFT of the stack of PSFs */
  result = psfstack_fft(psfn3, psfn2, psfn1, psfdata); /* Note this is an in-place transform! */
  if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(isverbose) fprintf(stdout, "--> Did FFT of PSF stack\n");
  psffftdata = psfdata; /* For naming convention's sake */
  if(stepwiseoutput) {
	result = write_minc_file("psffft.mnc", n3, n2, n1, psffftdata, COMPLEX_AS_COMPLEX);
	if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  }
  if(isverbose) fprintf(stdout, "--> Wrote file\n");

  /* Build the frequency-distance filter */
  result = build_fdp(psfn3,psfn2,psfn1, psffftdata, n3,n2,n1, fdpfilterdata, INVERSE_FILTER);
  if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(isverbose) fprintf(stdout, "--> Build FDP filter.\n");
  if(stepwiseoutput) {
	result = write_minc_file("fdpfilter.mnc", n3, n2, n1, fdpfilterdata, COMPLEX_AS_COMPLEX);
	if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  }
  if(isverbose) fprintf(stdout, "--> Wrote file\n");
  
  /* Build the noise dampening filter */
  result = build_noise_dampening_filter(n3,n2,n1,sgfftdata,noisedampendata);
  if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(isverbose) fprintf(stdout, "--> Built noise dampening filter.\n");
  if(stepwiseoutput) {
	result = write_minc_file("noisedampen.mnc", n3, n2, n1, noisedampendata, REAL_AS_REAL);
	if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  }
  if(isverbose) fprintf(stdout, "--> Wrote file\n");

  /* Build the rolloff filter if needed */

  /* Multiply the matrices to get final result */
  result = multiply_complex_by_real(nelems, sgfftdata, noisedampendata);
  if(!result) { free(sgfftdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(isverbose) fprintf(stdout, "--> Multiplied sg by noise dampen filter\n");
  if(stepwiseoutput) {
	result = write_minc_file("sgdampened.mnc", n3, n2, n1, sgfftdata, COMPLEX_AS_COMPLEX);
	if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  }
  if(isverbose) fprintf(stdout, "--> Wrote file\n");

  result = multiply_complex_by_complex(nelems, sgfftdata, fdpfilterdata);
  if(!result) { free(sgfftdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(isverbose) fprintf(stdout, "--> Completed matrix multiplication.\n");
  if(stepwiseoutput) {
	result = write_minc_file("filteredsg.mnc", n3, n2, n1, sgfftdata, COMPLEX_AS_COMPLEX);
	if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  }
  if(isverbose) fprintf(stdout, "--> Wrote file\n");


  /* Multiply by a Gaussian to deemphasize high-frequency */
  result = multiply_complex_by_gaussian(n3,n2,n1,sgfftdata,0.5);
  if(!result) { free(sgfftdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(stepwiseoutput) {
	result = write_minc_file("gauss_filtsg.mnc", n3, n2, n1, sgfftdata, COMPLEX_AS_COMPLEX);
	if(!result) { free(sgdata); free(psfdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  }
  if(isverbose) fprintf(stdout, "--> Multiplied by Gaussian roll off\n");

  /* Do the inverse transform back */
  result = sg_ccfft3d(n3,n2,n1,sgfftdata,INVERSE_FFT); /* Note this is an in-place transform! */
  if(!result) { free(sgfftdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(isverbose) fprintf(stdout, "--> Completed inverse 3DFFT\n");

  /* Write the output to file */
  result = write_minc_file(finalfile, n3, n2, n1, sgdata, COMPLEX_AS_REAL);
  if(!result) { free(sgfftdata); free(fdpfilterdata); free(noisedampendata); return 0; }
  if(isverbose) fprintf(stdout, "--> Wrote output to file.\n");

  /* Let's free up the memory we allocated */
  free(fdpfilterdata); free(noisedampendata); 
  free(sgdata); free(sgfftdata); 
  free(psfdata); free(psffftdata); /* Note we expect sgfft and psffft to be NULL */

  if(isverbose) fprintf(stdout,"Returning from main\n");
  return 1;

}

