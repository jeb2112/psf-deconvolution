
/******************************************
 *
 * mincsubvol
 *
 * Cuts down the minc file size by taking only
 * a subvolume of the entire data
 *
 * Johnathon Walls, 2006
 *
 ******************************************/


#include "fdr.h"

int isverbose=1;

int main(int argc, char *argv[]) {

  int result;
  long n3,n2,n1;
  long m3,m2,m1;
  int s0, s1, s2, c0, c1, c2;
  unsigned long start[3],count[3];
  float *data1 = NULL;
  
  if(argc!=5) {
	fprintf(stderr, "What are you doing?  I need a start, a count, an input file and an output file!\n");
	return 1;
  }

  if(! sscanf(argv[1], "%d,%d,%d",&s0,&s1,&s2)) {
    fprintf(stderr, "Start is improperly formatted!  (%s%)\n", argv[1]);
    return 0;
  }

  if(! sscanf(argv[2],"%d,%d,%d",&c0,&c1,&c2)) {
    fprintf(stderr, "Count is improperly formatted!  (%s%)\n", argv[2]);
    return 0;
  }

  fprintf(stdout, "Ok!  %s=%d-%d-%d   %s=%d-%d-%d\n", argv[1], s0,s1,s2,argv[2],c0,c1,c2);

  start[0] = s0;
  start[1] = s1;
  start[2] = s2;

  count[0] = c0;
  count[1] = c1;
  count[2] = c2;

  
  result = get_minc_subvol(argv[3], start, count, &n3, &n2, &n1, &data1, REAL_AS_REAL);
  if(!result) { free(data1); return 1; }
  if(isverbose) fprintf(stdout, "Loaded subvolume from minc file 1\n");

  result = write_minc_file(argv[4], count[0], count[1], count[2], data1, REAL_AS_REAL);
  if(isverbose) fprintf(stdout, "Wrote minc file\n");

  free(data1);
  
  return 0;

}
