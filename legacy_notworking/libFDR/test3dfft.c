
/******************************************
 *
 * forward3dfft
 *
 * Takes the magnitude of a complex data type minc2 file
 * and outputs it as a real minc2 file
 *
 * Johnathon Walls, 2005
 *
 ******************************************/


#include "fdr.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1;
  fdr_complex *data;
  
  if(argc!=3) {
	fprintf(stderr, "What are you doing?  I need an input file and an output file!\n");
	return 1;
  }

  result = open_minc_file(argv[1], &n3, &n2, &n1, &data, REAL_AS_COMPLEX);
  if(!result) { free(data); return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");

  /*
  data2 = calloc(n3*n2*n1,sizeof(fdr_complex));
  if(data2 == NULL) {
	fprintf(stderr, "Could not allocate memory for data2\n");
	return 1;
  }
  */
  
  result = sg_ccfft3d(n3,n2,n1,data,FORWARD_FFT);
  if(!result) return 0;
  if(isverbose) fprintf(stdout, "Calculated 3d FFT\n");
  
  result = sg_ccfft3d(n3,n2,n1,data,INVERSE_FFT);
  if(!result) return 0;

  result = write_minc_file(argv[2], n3, n2, n1, data, COMPLEX_AS_COMPLEX);
  if(isverbose) fprintf(stdout, "Wrote minc file\n");

  free(data);
  
  return 0;

}
