
/******************************************
 *
 * antialiasdata
 *
 * Takes an OPT data set and anti-aliases
 * along the view direction
 *
 * Johnathon Walls, 2006
 *
 ******************************************/

#define CUTOFFVALUE 550

#include "fdr2.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1, cutoff;
  mihandle_t filein, fileout;
  
  if(argc!=4) {
    fprintf(stderr, "What are you doing?  I need an input file and an output file and a cutoff value!\n");
    return 1;
  }

  result = open_minc_file_read(argv[1], &filein, COMPLEX);
  // result = open_minc_file(argv[1], &n3, &n2, &n1, &data, REAL_AS_COMPLEX);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Opened minc file 1 for reading\n");
  
  // Get dimensions from first file to duplicate for second file
  result = get_minc_dimensions_from_handle(filein, &n3, &n2, &n1);
  if(!result) return 1;
  
  // this should also check for clobber
  result = open_minc_file_write(argv[2], &fileout, n3*3, n2, n1, COMPLEX);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Opened minc file 2 for writing\n");
  
  cutoff = atof(argv[3]);
  result = slicewise_antialias_data(filein, fileout, cutoff);
  if(!result) return 1;
  if(isverbose) fprintf(stdout, "==>Completed antialias\n");
  
  miclose_volume(filein);
  miclose_volume(fileout);
  if(isverbose) fprintf(stdout, "==>Closed all minc files\n");

  return 0;


/*
  int result;
  long n3,n2,n1;
  long cutoff;
  fdr_complex *fftdata, *outdata;
  
  if(argc!=4) {
	fprintf(stderr, "What are you doing?  I need an input file, a cutoff value, and an output file!\n");
	return 1;
  }

  result = open_minc_file(argv[1], &n3, &n2, &n1, &fftdata, COMPLEX_AS_COMPLEX);
  if(!result) { free(fftdata); return 1; }
  if(isverbose) fprintf(stdout, "Loaded minc file\n");
  
  cutoff = atof(argv[2]);
  outdata = calloc(n3*3*n2*n1, sizeof(fdr_complex));
  if(outdata == NULL) {
	fprintf(stderr, "Could not allocate memory for outdata!\n");
	free(fftdata); free(outdata);
	return 1;
  }
	
  result = build_antialiased_data(n3,n2,n1, cutoff, fftdata, outdata);
  if(!result) { free(fftdata); free(outdata); return 1; }
  if(isverbose) fprintf(stdout, "Anti-aliased data.\n");
  
  // print_data(n3*3,n2,n1,outdata);

  result = write_minc_file(argv[3], n3*3, n2, n1, outdata, COMPLEX_AS_COMPLEX);
  if(isverbose) fprintf(stdout, "Wrote minc file\n");

  free(fftdata); free(outdata);
  
  return 0;*/

}
