
/******************************************
 *
 * buildbwfilter
 *
 * Calculates the rolloff filter for FDR
 *
 * Johnathon Walls, 2005
 *
 ******************************************/


#include "fdr2.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1;
  double bandlimit;
  mihandle_t templatefile, fileout;
  
  if(argc!=4) {
    fprintf(stderr, "What are you doing?  I need a template file, a bandlimit, and an output file!\n");
    return 1;
  }

  bandlimit = atof(argv[2]);
  
  // fprintf(stdout, "Proceeding with bandlimit of %2.2f\n", bandlimit);
  
  result = open_minc_file_read(argv[1], &templatefile, FLOAT);
  if(!result) return 1;
  if(isverbose) fprintf(stdout,"==>Opened minc template file.\n");

  result = get_minc_dimensions_from_handle(templatefile, &n3, &n2, &n1);
  if(!result) return 1;
  if(isverbose) fprintf(stdout,"==>Retrieved dimensions from template.\n");
   
  result = open_minc_file_write(argv[3], &fileout, n3, n2, n1, FLOAT);
  if(!result) return 1;   
  if(isverbose) fprintf(stdout, "==>Opened file for output\n");

  // fprintf(stdout, "Bandlimit is %2.2f\n", bandlimit);
  
  result = slicewise_build_bw_filter(fileout, bandlimit);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Built bandlimit filter.\n");
  
  return 0;

}
