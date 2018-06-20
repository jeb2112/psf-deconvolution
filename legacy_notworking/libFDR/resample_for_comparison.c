
/******************************************
 *
 * resample_for_comparison
 *
 * Loads in the 3d fft of a sinogram and converts
 * it into the psf-stack space
 *
 * Johnathon Walls, 2006
 *
 ******************************************/


#include "fdr.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1;
  long psfn3, psfn2, psfn1;
  fdr_complex *sgfftdata = NULL;
  fdr_complex *psfcomparedata = NULL;
  
  if(argc!=4) {
	fprintf(stderr, "Syntax: buildcomparefdp <psf_template> <input> <output>\n");
	return 1;
  }


  result = get_minc_dimensions(argv[1], &psfn3, &psfn2, &psfn1);
  if(!result) return 0;

  result = open_minc_file(argv[2], &n3, &n2, &n1, &sgfftdata, COMPLEX_AS_COMPLEX);
  if(!result) { free(sgfftdata); return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");

  if(n2 != psfn2 || n1 != psfn1) {
	fprintf(stderr, "X-Y dimensions do not match up (%d,%d) vs (%d,%d)\n", n2, n1, psfn2, psfn1);
	free(sgfftdata);
	return 0;
  }
  if(isverbose) fprintf(stdout, "Dimensions match up, that's good.\n");

  psfcomparedata = (fdr_complex *)malloc(psfn3*psfn2*psfn1*sizeof(float));
  if(psfcomparedata == NULL) {
	fprintf(stderr, "Could not allocate memory for psfcomparedata!\n");
	free(sgfftdata); free(psfcomparedata);
  }
  if(isverbose) fprintf(stdout, "Allocated memory for fdpfilterdata.\n");

  result = resample_for_comparison(psfn3,psfn2,psfn1, psfcomparedata, 
                                   n3,n2,n1, sgfftdata);
  if(!result) { free(sgfftdata); free(psfcomparedata); return 0; }
  if(isverbose) fprintf(stdout, "Built FDP comparison filter.\n");

  result = write_minc_file(argv[3], psfn3, psfn2, psfn1, psfcomparedata, REAL_AS_REAL);
  if(isverbose) fprintf(stdout, "Wrote minc file\n");

  free(psfcomparedata); free(sgfftdata);
  
  return 0;

}
