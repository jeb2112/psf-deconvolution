
/******************************************
 *
 * psfstackfft
 *
 *
 * Johnathon Walls, 2006
 *
 ******************************************/


#include "fdr2.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1;
  mihandle_t psffile, normfile;

  if(argc!=3) {
    fprintf(stderr, "What are you doing?  I need an input file and an output file!\n");
    return 1;
  }

  result = open_minc_file_read(argv[1], &psffile, FLOAT);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Opened minc file for reading\n");

  result = get_minc_dimensions_from_handle(psffile, &n3, &n2, &n1);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Opened minc file for writing\n");
  
  result = open_minc_file_write(argv[2], &normfile, n3, n2, n1, FLOAT);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Opened minc file for reading\n");

  result = slicewise_normalize_psfstack(psffile, normfile);
  if(!result) { return 1; }
  if(isverbose) fprintf(stdout, "==>Completed normalization of PSF stacks.\n");
  
  return 0;

}
