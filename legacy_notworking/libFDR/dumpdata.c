
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
  long n3,n2,n1,i;
  fdr_complex *data;
  
  if(argc!=2) {
	fprintf(stderr, "What are you doing?  I need an input file!\n");
	return 1;
  }

  result = open_minc_file(argv[1], &n3, &n2, &n1, &data, REAL_AS_REAL);
  if(!result) { free(data); return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");

  for(i=0;i<n3*n2*n1;i++) {
	fprintf(stdout, "Value at index %i is %f %f\n", i, data[i].re, data[i].im);
  }

  free(data);
  
  return 0;

}
