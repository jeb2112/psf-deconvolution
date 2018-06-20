
/******************************************
 *
 * forward3dfft
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
  mihandle_t filein, fileout;
  
  if(argc!=3) {
	fprintf(stderr, "What are you doing?  I need an input file and an output file!\n");
	return 1;
  }

	result = open_minc_file_read(argv[1], &filein, FLOAT);
  // result = open_minc_file(argv[1], &n3, &n2, &n1, &data, REAL_AS_COMPLEX);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Opened minc file 1 for reading\n");
  
  // Get dimensions from first file to duplicate for second file
  
  result = get_minc_dimensions_from_handle(filein, &n3, &n2, &n1);
	if(!result) return 1;
  
  // this should also check for clobber
  result = open_minc_file_write(argv[2], &fileout, n3, n2, n1, COMPLEX);
	if(!result) return 1;
	if(isverbose) fprintf(stdout, "==>Opened minc file 2 for writing\n");
	
	result = full_2dfft(filein, fileout, 0);
  // result = sg_ccfft3d(n3,n2,n1,data,FORWARD_FFT);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Calculated 3d FFT\n");
  
  // result = write_minc_file(argv[2], n3, n2, n1, data, COMPLEX_AS_COMPLEX);
  miclose_volume(filein);
  miclose_volume(fileout);
  if(isverbose) fprintf(stdout, "==>Closed all minc files\n");

  return 0;

}
