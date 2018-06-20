
/******************************************
 *
 * buildinversefdp
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
  long n3,n2,n1;
  double limit;
  mihandle_t infile, outfile;

  if(argc!=4) {
    fprintf(stderr, "What are you doing?  I need an input file, a limit value, and an output file!\n");
    return 1;
  }

  limit = atof(argv[2]);

  result = open_minc_file_read(argv[1], &infile, COMPLEX);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Loaded minc file\n");
  
  result = get_minc_dimensions_from_handle(infile, &n3, &n2, &n1);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Retrieved dimensions\n");

  result = open_minc_file_write(argv[3], &outfile, n3, n2, n1, COMPLEX);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Created minc file 2\n");

	result = slicewise_limit_recovery_filter(infile, outfile, limit);
	if(!result) { return 1; }
 	if(isverbose) fprintf(stdout, "==>Limited FDR filter.\n");

  
  return 0;

}
