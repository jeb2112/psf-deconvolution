
/******************************************
 *
 * buildcomparefdp
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
  long psfn3, psfn2, psfn1;
  float maxslope, bw;
  fdr_complex *fdpfilterdata = NULL;
  fdr_complex *psfdata = NULL;
  
  if(argc!=6) {
	fprintf(stderr, "Syntax: buildcomparefdp <template> <maxslope> <bandwidth> <input> <output>\n");
	return 1;
  }

  maxslope = atof(argv[2]);
  bw = atof(argv[3]);


  result = get_minc_dimensions(argv[1], &n3, &n2, &n1);
  if(!result) return 0;

  result = open_minc_file(argv[4], &psfn3, &psfn2, &psfn1, &psfdata, COMPLEX_AS_COMPLEX);
  if(!result) { free(psfdata); return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");

  if(n2 != psfn2 || n1 != psfn1) {
	fprintf(stderr, "X-Y dimensions do not match up (%d,%d) vs (%d,%d)\n", n2, n1, psfn2, psfn1);
	free(psfdata);
	return 0;
  }
  if(isverbose) fprintf(stdout, "Dimensions match up, that's good.\n");

  fdpfilterdata = (fdr_complex *)malloc(n3*n2*n1*sizeof(fdr_complex));
  if(fdpfilterdata == NULL) {
	fprintf(stderr, "Could not allocate memory for fdpfilterdata!\n");
	free(psfdata); free(fdpfilterdata);
  }
  if(isverbose) fprintf(stdout, "Allocated memory for fdpfilterdata.\n");

  result = build_fdp(psfn3,psfn2,psfn1, psfdata, n3,n2,n1, fdpfilterdata, maxslope, bw, COMPARE_FILTER);
  if(!result) { free(psfdata); free(fdpfilterdata); return 0; }
  if(isverbose) fprintf(stdout, "Built FDP comparison filter.\n");

  result = write_minc_file(argv[5], n3, n2, n1, fdpfilterdata, COMPLEX_AS_COMPLEX);
  if(isverbose) fprintf(stdout, "Wrote minc file\n");

  free(psfdata); free(fdpfilterdata);
  
  return 0;

}
