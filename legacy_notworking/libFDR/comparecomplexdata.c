
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


#include "fdr2.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  mihandle_t file1, file2;

  if(argc!=3) {
	fprintf(stderr, "What are you doing?  I need file1 and file2!\n");
	return 1;
  }
  
  fprintf(stdout, "Beginning comparison ...\n");
  result = open_minc_file_read(argv[1], &file1, COMPLEX);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file 1\n");

  result = open_minc_file_read(argv[2], &file2, COMPLEX);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");

  result = slicewise_compare_complex_data(file1, file2);
  if(!result) { return 1; }

  fprintf(stdout, "Comparison complete!\n");
  
  return 0;

}
