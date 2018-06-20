
/* #define SGI 0 */
#define LINUX 1
/* #define OMP 0 */

#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include "minc2.h"
/*
#ifdef OMP
#    include "omp.h"
#endif
*/
#define COMPLEX_AS_COMPLEX  3
#define COMPLEX_AS_REAL     2
#define REAL_AS_COMPLEX     1
#define REAL_AS_REAL        0
#define FORWARD_FFT         1
#define INVERSE_FFT        -1
#define COMPARE_FILTER      0
#define INVERSE_FILTER      1

/* SGI */
/*
#include "scsl_fft"
typedef scsl_complex fdr_complex;
*/
/* Linux */
#include "complex.h"
#include "fftw3.h"
typedef fftw_complex fdr_complex;

/* These functions are used for low-memory systems */
int print_data(long n3, long n2, long n1, fdr_complex* data);
float sgsum(int elems, float *data) ;
float mymax(int elems, float *data);
float mymin(int elems, float *data);
float absmax(int elems,fdr_complex *data);
float complexmax(int elems,fdr_complex *data);
float complexmin(int elems,fdr_complex *data);
int correct_fft_phaseshift(long Z, long Y, long X, float *data);
int shift_complex_to_float(long elems, float *data);
int shift_float_to_complex(long elems, fdr_complex *data);
int calculate_abs(long nelems, fdr_complex *data);
int calculate_phase(long nelems, fdr_complex *data);
int open_minc_file(char *filein, long *nz, long *ny, long *nx, void **data, int flags);
int get_minc_subvol(char *filein, unsigned long *start, unsigned long *count, long *nz, long *ny, long *nx, void **data, int flags);
int get_minc_dimensions(char *filein, long *nz, long *ny, long *nx);
int write_minc_file(char *fileout, long nz, long ny, long nx, void *data, int flags);
int cfftshift(long Z,long Y,long X, fdr_complex* data);
int c2dfftshift(long Y,long X, fdr_complex* data);
int fftshift(long Z,long Y,long X, float* data);
int sg_ccfft3d(long n3, long n2, long n1, fdr_complex *data, int direction);
int psfstack_fft(long n3, long n2, long n1, fdr_complex *data);
int normalize_psfstack(long n3, long n2, long n1, float *data);
int normalize_complexpsfstack(long n3, long n2, long n1, fdr_complex *data);
int limit_recovery_filter(long nelems, float limit, fdr_complex *filterdata);
int compare_complex_data(long n3, long n2, long n1, fdr_complex *data1, fdr_complex *data2);
int compare_real_data(long n3, long n2, long n1, float *data1, float *data2);
int build_rolloff_filter(long n3, long n2, long n1, float weight, float maxslope, float *outdata);
int build_wiener_filter(long n3, long n2, long n1, long newelems, fdr_complex *fftdata, float *outdata);
int build_noise_dampening_filter(long n3, long n2, long n1, fdr_complex *fftdata, float *outdata);
int build_bw_filter(long n3, long n2, long n1, float bandlimit, float *bwfilterdata);
int resample_for_comparison(long psfn3, long psfn2, long psfn1, float *psfcomparedata,
                               long    n3, long    n2, long    n1, fdr_complex *sgfftdata);
int build_fdp(long psfn3, long psfn2, long psfn1, fdr_complex *psffftdata, 
			  long    n3, long    n2, long    n1, fdr_complex *fdpfilterdata, 
			  float mxsl, float bandlimit, int flags);
int multiply_complex_by_complex(long nelems, fdr_complex *data1, fdr_complex *data2);
int multiply_complex_by_real(long nelems, fdr_complex *data1, float *data2);
int multiply_complex_by_gaussian(long n3, long n2, long n1, fdr_complex *data, double weight);
int multiply_real_by_gaussian(long n3, long n2, long n1, float *data, double weight);
int build_antialiased_data(long n3, long n2, long n1, long cutoff, fdr_complex *fftdata, fdr_complex *outdata);

/* These functions are used only for large memory systems */
/*
int print_data(long n3, long n2, long n1, fdr_complex* data);
float sgsum(int elems, float *data) ;
float mymax(int elems, float *data);
float mymin(int elems, float *data);
float absmax(int elems,fdr_complex *data);
float complexmax(int elems,fdr_complex *data);
float complexmin(int elems,fdr_complex *data);
int correct_fft_phaseshift(long Z, long Y, long X, float *data);
int shift_complex_to_float(long elems, float *data);
int shift_float_to_complex(long elems, fdr_complex *data);
int calculate_abs(long nelems, fdr_complex *data);
int calculate_phase(long nelems, fdr_complex *data);
int open_minc_file(char *filein, long *nz, long *ny, long *nx, void **data, int flags);
int get_minc_subvol(char *filein, unsigned long *start, unsigned long *count, long *nz, long *ny, long *nx, void **data, int flags);
int get_minc_dimensions(char *filein, long *nz, long *ny, long *nx);
int write_minc_file(char *fileout, long nz, long ny, long nx, void *data, int flags);
int cfftshift(long Z,long Y,long X, fdr_complex* data);
int c2dfftshift(long Y,long X, fdr_complex* data);
int fftshift(long Z,long Y,long X, float* data);
int sg_ccfft3d(long n3, long n2, long n1, fdr_complex *data, int direction);
int psfstack_fft(long n3, long n2, long n1, fdr_complex *data);
int normalize_psfstack(long n3, long n2, long n1, float *data);
int normalize_complexpsfstack(long n3, long n2, long n1, fdr_complex *data);
int limit_recovery_filter(long nelems, float limit, fdr_complex *filterdata);
int compare_complex_data(long n3, long n2, long n1, fdr_complex *data1, fdr_complex *data2);
int compare_real_data(long n3, long n2, long n1, float *data1, float *data2);
int build_rolloff_filter(long n3, long n2, long n1, float weight, float maxslope, float *outdata);
int build_wiener_filter(long n3, long n2, long n1, long newelems, fdr_complex *fftdata, float *outdata);
int build_noise_dampening_filter(long n3, long n2, long n1, fdr_complex *fftdata, float *outdata);
int build_bw_filter(long n3, long n2, long n1, float bandlimit, float *bwfilterdata);
int resample_for_comparison(long psfn3, long psfn2, long psfn1, float *psfcomparedata,
                               long    n3, long    n2, long    n1, fdr_complex *sgfftdata);
int build_fdp(long psfn3, long psfn2, long psfn1, fdr_complex *psffftdata, 
			  long    n3, long    n2, long    n1, fdr_complex *fdpfilterdata, 
			  float mxsl, float bandlimit, int flags);
int multiply_complex_by_complex(long nelems, fdr_complex *data1, fdr_complex *data2);
int multiply_complex_by_real(long nelems, fdr_complex *data1, float *data2);
int multiply_complex_by_gaussian(long n3, long n2, long n1, fdr_complex *data, double weight);
int multiply_real_by_gaussian(long n3, long n2, long n1, float *data, double weight);
int build_antialiased_data(long n3, long n2, long n1, long cutoff, fdr_complex *fftdata, fdr_complex *outdata);
*/