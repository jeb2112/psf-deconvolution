
/******************************************
 *
 * comparecomplexdata
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
  float *data1 = NULL;
  float *data2 = NULL;

  if(argc!=3) {
	fprintf(stderr, "What are you doing?  I need a lower valued file and a higher valued file!\n");
	return 1;
  }
  
  fprintf(stdout, "Beginning comparison ...\n");
  result = open_minc_file(argv[1], &n3, &n2, &n1, &data1, REAL_AS_REAL);
  if(!result) { free(data1); return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");

  result = open_minc_file(argv[2], &n3, &n2, &n1, &data2, REAL_AS_REAL);
  if(!result) { free(data1); free(data2); return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");

  result = compare_real_data(n3,n2,n1, data1, data2);
  if(!result) { free(data1); free(data2); return 1; }

  fprintf(stdout, "Comparison complete!\n");

  free(data1); free(data2);
  
  return 0;

}
