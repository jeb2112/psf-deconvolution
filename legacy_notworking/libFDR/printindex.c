
/******************************************
 *
 * printrealvalue
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
  long index;
  float *data = NULL;

  if(argc!=3) {
	fprintf(stderr, "What are you doing?  I need a lower valued file and an integer!\n");
	return 1;
  }
  
  fprintf(stdout, "Beginning comparison ...\n");
  result = open_minc_file(argv[1], &n3, &n2, &n1, &data, REAL_AS_REAL);
  if(!result) { free(data); return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");

  index = atol(argv[2]);
  if(index > n3*n2*n1-1) {
	fprintf(stderr, "Index %d is larger than max possible index %d\n", index, n3*n2*n1-1);
	free(data);
  }

  fprintf(stdout, "The value at index %d is %f

  free(data);
  
  return 0;

}
