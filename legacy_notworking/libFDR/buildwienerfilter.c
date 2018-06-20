
/******************************************
 *
 * buildnoisedampenfilter
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
  mihandle_t filein, fileout;
  long n3,n2,n1;
  
  if(argc!=3) {
	fprintf(stderr, "What are you doing?  I need an input file and an output file!\n");
	return 1;
  }

  result = open_minc_file_read(argv[1], &filein, FLOAT);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Opened minc file 1 for reading\n");

  result = get_minc_dimensions_from_handle(filein, &n3, &n2, &n1);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Got dimensions from minc file 1\n");
  
  result = open_minc_file_write(argv[2], &fileout, n3, n2, n1, FLOAT);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Opened minc file 2 for writing\n");

  result = slicewise_build_wiener_filter(filein, fileout, 500);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Built wiener filter.\n");
  
  return 0;

}
