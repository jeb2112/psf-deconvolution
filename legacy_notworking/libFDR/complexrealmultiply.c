
/******************************************
 *
 * complexmultiply
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

  mihandle_t file1, file2, fileout;
  
  if(argc!=4) {
    fprintf(stderr, "What are you doing?  I need two input files and an output file!\n");
    return 1;
  }

  result = open_minc_file_read(argv[1], &file1, COMPLEX);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==> Opened file 1\n");
  
  result = open_minc_file_read(argv[2], &file2, FLOAT);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==> Opened file 2\n");
    
  // Get dimensions from first file to duplicate for second file
  
  result = get_minc_dimensions_from_handle(file1, &n3, &n2, &n1);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==> Obtained dimensions\n");
    
  open_minc_file_write(argv[3], &fileout, n3, n2, n1, COMPLEX);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==> Opened file out\n");
    
  result = slicewise_multiply_complex_by_real(file1, file2, fileout, 0);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "Multiplied data sets\n");
  
  miclose_volume(file1); miclose_volume(file2); miclose_volume(fileout);
  
  return 0;

}
