
/******************************************
 *
 * rfftshift
 *
 * Shifts the fft data
 *
 * Johnathon Walls, 2005
 *
 ******************************************/


#include "fdr2.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1;

  mihandle_t filein, fileout;
  
  if(argc!=3) {
    fprintf(stderr, "What are you doing?  I need an input file and an output file!\n");
    return 1;
  }

  result = open_minc_file_read(argv[1], &filein, FLOAT);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==> Opened file in\n");
      
  // Get dimensions from first file to duplicate for second file
  
  result = get_minc_dimensions_from_handle(filein, &n3, &n2, &n1);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==> Obtained dimensions\n");
    
  open_minc_file_write(argv[2], &fileout, n3, n2, n1, FLOAT);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==> Opened file out\n");
    
  result = slicewise_rfftshift(filein, fileout);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "Completed FFT shift\n");
  
  miclose_volume(filein);  miclose_volume(fileout);
  
  return 0;

}
