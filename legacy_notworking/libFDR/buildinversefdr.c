
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
  long psfn3, psfn2, psfn1;
  float limit, slope, bw;
  mihandle_t template_file, filter_file, psf_file;
  
  limit = -1;

  if(argc!=4) {
  	fprintf(stderr, "Syntax: buildinversefdr <template> <psffile> <outfile>\n");
  	return 1;
  }

  result = open_minc_file_read(argv[1], &template_file, COMPLEX);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Opened template file\n");
  
  result = get_minc_dimensions_from_handle(template_file, &n3, &n2, &n1);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Retrieved dimensions\n");

  result = open_minc_file_read(argv[2], &psf_file, COMPLEX);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Opened psf file\n");
  
  result = get_minc_dimensions_from_handle(psf_file, &psfn3, &psfn2, &psfn1);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Retrieved dimensions from psf file\n");
  

  if(n2 != psfn2 || n1 != psfn1) {
  	fprintf(stderr, "X-Y dimensions do not match up (%d,%d) vs (%d,%d)\n", n2, n1, psfn2, psfn1);
  	return 1;
  }
  if(isverbose) fprintf(stdout, "==>Dimensions match up, that's good.\n");

  result = open_minc_file_write(argv[3], &filter_file, n3, n2, n1, COMPLEX);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Opened filter file for writing.\n");

  result = slicewise_build_inverse_fdr_filter(psf_file, filter_file, INVERSE_FILTER);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Built Inverse FDR filter.\n");
 
  return 0;

}
