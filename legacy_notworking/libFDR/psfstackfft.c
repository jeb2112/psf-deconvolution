
/******************************************
 *
 * psfstackfft
 *
 *
 * Johnathon Walls, 2006
 *
 ******************************************/


#include "fdr2.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1;
  mihandle_t psf_file, psffft_file;
  

  if(argc!=3) {
	  fprintf(stderr, "What are you doing?  I need an input file and an output file!\n");
	  return 1;
  }

  result = open_minc_file_read(argv[1], &psf_file, FLOAT);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Opened PSF file\n");
  
  result = get_minc_dimensions_from_handle(psf_file, &n3, &n2, &n1);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Retrieved PSF dimensions\n");

  result = open_minc_file_write(argv[2], &psffft_file, n3, n2, n1, COMPLEX);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Opened PSF FFT file for writing\n");

  result = slicewise_psfstack_fft(psf_file, psffft_file);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Completed FFT of PSF stacks.\n");

  return 0;

}
