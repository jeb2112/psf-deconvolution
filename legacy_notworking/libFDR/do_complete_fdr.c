
/******************************************
 *
 * do_complete_fdr
 *
 * Does everything, from start to finish
 *
 * Johnathon Walls, 2006
 *
 ******************************************/

#include "fdr.h"

#define WIENER_ELEMENTS 500

int isverbose=1;

int clean_filter(long nelems, float *data) {

  int i;

  for(i=0;i<nelems;i++) {
    data[i] = 0;
  }

  return 1;

}

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1, psfn3, psfn2, psfn1;
  float slope, bw, rolloff, limit;
  float *buffer1;
  fdr_complex *data, *buffer2, *psf;
  
  /* Arguments are formatted as: 
   * do_complete_fdr SLOPE BANDLIMIT ROLLOFF FDRLIMIT data_file psf_file out_file
   */

  if(argc!=8) {
	fprintf(stderr, "Syntax is: do_complete_fdr SLOPE BANDLIMIT ROLLOFF FDRLIMIT <data_file> <psf_file> <out_file>\n");
	return 1;
  }

  slope = atof(argv[1]);
  bw = atof(argv[2]);
  rolloff = atof(argv[3]);
  limit = atof(argv[4]);

  /* 1 - Open up the data we need */
  result = open_minc_file(argv[5], &n3, &n2, &n1, &data, REAL_AS_COMPLEX);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Loaded minc file %s\n", argv[5]);
  
  /* 2 - Open up PSF file */
  result = open_minc_file(argv[6], &psfn3, &psfn2, &psfn1, &psf, REAL_AS_COMPLEX);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Loaded PSF minc file %s\n", argv[6]);

  /* 3 - Check dimensions */
  if(n2 != psfn2 || n1 != psfn1) {
        fprintf(stderr, "X-Y dimensions do not match up (%d,%d) vs (%d,%d)\n", n2, n1, psfn2, psfn1);
        free(data); free(data); free(buffer1); free(buffer2);
        return 0;
  }
  if(isverbose) fprintf(stdout, "==> Dimensions match up, that's good.\n");

  /* 4 - Allocate additional storage space */
  buffer1 = (float *)malloc(n3*n2*n1*sizeof(float));
  if(buffer1 == NULL) {
        fprintf(stderr, "Could not allocate memory for float buffer!\n");
        free(data); free(psf); free(buffer1); free(buffer2);
  }
  if(isverbose) fprintf(stdout, "==> Allocated memory for buffer1.\n");

  /* 5 - Build wiener filter */
  result = build_wiener_filter(n3,n2,n1, WIENER_ELEMENTS, data, buffer1);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Built wiener filter.\n");

  /* 6 - Multiply data by wiener filter */  
  result = multiply_complex_by_real(n3*n2*n1,data, buffer1);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Multiplied data sets\n");

  /* 7 - Build bandwidth filter */
  result = build_bw_filter(n3,n2,n1, bw, buffer1);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Built bw filter.\n");

  /* 8 - multiply data by bandwidth fitler */
  clean_filter(n3*n2*n1, buffer1);
  result = multiply_complex_by_real(n3*n2*n1,data, buffer1);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Multiplied data sets\n");

  /* 9 - build rolloff filter */
  result = build_rolloff_filter(n3,n2,n1, rolloff, slope, buffer1);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Built rolloff filter.\n");

  /* 10 - multiply data by rolloff filter */
  clean_filter(n3*n2*n1, buffer1);
  result = multiply_complex_by_real(n3*n2*n1,data, buffer1);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Multiplied data sets\n");

  /* 11 - free float buffer */
  free(buffer1);

  /* 12 - Allocate complex storage space */
  buffer2 = (fdr_complex *)malloc(n3*n2*n1*sizeof(fdr_complex));
  if(buffer2 == NULL) {
        fprintf(stderr, "Could not allocate memory for fdr_complex buffer!\n");
        free(data); free(psf); free(buffer1); free(buffer2);
  }
  if(isverbose) fprintf(stdout, "==> Allocated memory for buffer2.\n");

  /* 13 - Normalize psf */
  result = normalize_complexpsfstack(psfn3, psfn2, psfn1, psf);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Normalized PSF\n");

  /* 14 - PSF FFT along distance axis */
  result = psfstack_fft(psfn3, psfn2, psfn1, (fdr_complex *)psf);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Completed FFT of PSF stacks.\n");

  /* 15 - Build inverse FDR filter */
  result = build_fdp(psfn3,psfn2,psfn1, psf, n3,n2,n1, buffer2, slope, bw, INVERSE_FILTER);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Built FDR filter.\n");
                                                                               
  /* 16 - Limit FDR */
  result = limit_recovery_filter(n3*n2*n1, limit, buffer2);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Limited FDR filter.\n");

  /* 17 - Forward FFT of views */
  result = sg_ccfft3d(n3,n2,n1,data,FORWARD_FFT);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Calculated 3d FFT\n");

  /* 18 - Multiply data fft by inverse FDR */
  result = multiply_complex_by_complex(n3*n2*n1,data,buffer2);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Multiplied data FFT by FDR\n");

  /* 19 - Free the complex data buffer - we don't need it anymore */
  free(buffer2);

  /* 20 - Calculate the magnitude of the end result! */
  result = calculate_abs(n3*n2*n1, data);
  if(!result) { free(data); free(psf); free(buffer1); free(buffer2); return 0; }
  if(isverbose) fprintf(stdout, "==> Multiplied data FFT by FDR\n");
  
  /* 21 - Write out the final result! */
  result = write_minc_file(argv[7], n3, n2, n1, data, REAL_AS_REAL);
  if(isverbose) fprintf(stdout, "==> Wrote minc file as %s\n", argv[7]);

  free(data); free(psf); free(buffer1); free(buffer2);
  
  return 0;

}
