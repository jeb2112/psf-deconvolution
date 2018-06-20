
/******************************************
 *
 * buildrollofffilter
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
  double weight;
  double maxslope;
  mihandle_t templatefile, fileout;
  
  if(argc!=5) {
	 fprintf(stderr, "What are you doing?  I need a template file, a weight, a maxslope, and an output file!\n");
	 return 1;
  }

  weight = atof(argv[2]);
  maxslope = atof(argv[3]);

  result = open_minc_file_read(argv[1], &templatefile, FLOAT);
  if(!result) return 1;
  if(isverbose) fprintf(stdout,"==>Opened minc template file.\n");

  result = get_minc_dimensions_from_handle(templatefile, &n3, &n2, &n1);
  if(!result) return 1;
  if(isverbose) fprintf(stdout,"==>Retrieved dimensions from template.\n");
   
  result = open_minc_file_write(argv[4], &fileout, n3, n2, n1, FLOAT);
  if(!result) return 1;   
  if(isverbose) fprintf(stdout, "==>Opened file for output\n");

  result = slicewise_build_rolloff_filter(fileout, weight, maxslope);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Built rolloff filter.\n");

  return 0;

}
