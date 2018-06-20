#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include "math.h"

typedef fftwf_complex fdrf_complex;
typedef fftw_complex fdr_complex;

int X=1000;
int Y=1000;
int Z=200;
double precision = 1e-8;

inline int myrand() {
  return rand()/(int)(((unsigned)RAND_MAX+1)/4095);
}

int twodtest()
{
	  //float sgdata[Z][Y][X];
  //complex sgfftdata[Z][Y][X/2];
  //float table[(Z + 256) + (2*Y + 256) + (2*X + 256)];
  //float work[Z + 4*X];
  fdr_complex* orig = NULL;
  fdr_complex* origtracker = NULL;
  fdr_complex* fft = NULL;
  fdr_complex* back = NULL;
  fftw_plan p1,p2;
  long x,y,z;
  long index, index2;
  long total = 0;
  double complex sum = 0;

  orig = (fdr_complex *) fftw_malloc( X * Y * sizeof(fdr_complex));
  fft  = (fdr_complex *) fftw_malloc( X * Y * sizeof(fdr_complex));
  back = (fdr_complex *) fftw_malloc( X * Y * sizeof(fdr_complex));

  /* Build the plan */
  // fftw_plan_with_nthreads(1);
  p1 = fftw_plan_dft_2d(X,Y, orig, fft, FFTW_FORWARD, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_2d(X,Y, fft, back, FFTW_BACKWARD, FFTW_ESTIMATE);

  printf("plan built\n");


	for(z=0;z<Z;z++)
	{
		for (y=0;y<Y;y++) {
	  	for (x=0;x<X;x++) {
			index = y*X + x;
		  	orig[index] = 1.0*myrand() + 0.5;
		  	sum += orig[index];
          //          printf("At %d we have %f %f (from %f)\n", index, orig[index][0], orig[index][1], 1.0*myrand() + 0.5);
	  	}
		}
		// printf("values set\n");
		
	  fftw_execute(p1);

  /*
  isys[0] = 1;
  zzfft3d(0, X, Y, Z, 1.0f, (complex *)sgdata, X, Y, (complex *)sgfftdata, X, Y,
		  table, work, isys);
  zzfft3d(1, X, Y, Z, 1.0f, (complex *)sgdata, X, Y, (complex *)sgfftdata, X, Y,
          table, work, isys);
  */
  	// printf("Forward FFT complete\n");


	  fftw_execute(p2);
  //  csfft3d(-1, X, Y, Z, 1.0f/(1.0*X*Y*Z),(complex *) sgfftdata, X2, Y, 
  //		  (float *) sgdata, X+2, Y,table, work, isys);

		// printf("Firsts are: %f %f %f %f\n", orig[0], back[0]);
   	for(y=0;y<Y;y++)
   	{
   		// printf("Y IS-----------------------> %d\n", y);
   		for(x=0;x<X;x++)
   		{
   			// printf("%d\n",x);
   			index = y*X + x;
   			back[index] /= Y*X;
   			if( (fabs( (orig[index]-back[index])/(orig[index])) > precision) && creal(orig[index])) 
   			{
					// printf("Values are %d : %f %f\n", (int)index, orig[index],back[index2]);
		  		total++;
   			}
   		}
  	}
  	// printf("Backward FFT complete!\n");
	}
/*
  //#pragma omp parallel for
  for (z=0;z<Z;z++) {
		for (y=0;y<Y;y++) {
	  	for (x=0;x<X;x++) {
				index = z*Y*X + y*X + x;
        back[index] /= X*Y*Z;
        // back[index] /= X*Y*Z;

				if( (fabs( (orig[index]-back[index])/(orig[index])) > precision) && orig[index]) {
					printf("Values are %d : %f %f\n", (int)index, orig[index],back[index]);
		  		total++;
				}
	  	}
		}
  }
*/
  printf("[Double] Total number of voxels less than precision (%f %%): %d\n", precision*100, total);
  printf("[Double] Sum is: %ld, making a mean of %f\n", (long)(sum), (long)(sum)/((float)X*Y*Z));
  printf("TESTING %f %f\n", fft[85]);
/*
  for (z=0;z<Z;z++) {
    for (y=0;y<Y;y++) {
      for (x=0;x<X;x++) {
  		index = z*Y*X + y*X + x;
        //        printf("Values are %d : %f %f\n", (int)index, orig[index][0],back[index][0]);
  	  }
  	}
  }
  */
  //printf("End expected at %d\n", Z*Y*X2P1-1);

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_free(orig); fftw_free(fft); fftw_free(back); fftw_free(origtracker);
  // fftw_cleanup_threads();


  return 1;
}

int doubletest() {

  //float sgdata[Z][Y][X];
  //complex sgfftdata[Z][Y][X/2];
  //float table[(Z + 256) + (2*Y + 256) + (2*X + 256)];
  //float work[Z + 4*X];
  fdr_complex* orig = NULL;
  fdr_complex* fft = NULL;
  fdr_complex* back = NULL;
  fftw_plan p1,p2;
  long x,y,z;
  long index;
  long total = 0;
  double complex sum = 0;

  orig = (fdr_complex *) fftw_malloc( X * Y * Z * sizeof(fdr_complex));
  fft  = (fdr_complex *) fftw_malloc( X * Y * Z * sizeof(fdr_complex));
  back = (fdr_complex *) fftw_malloc( X * Y * Z * sizeof(fdr_complex));

  /* Build the plan */
  // fftw_plan_with_nthreads(1);
  p1 = fftw_plan_dft_3d(X,Y,Z, orig, fft, FFTW_FORWARD, FFTW_ESTIMATE);
  p2 = fftw_plan_dft_3d(X,Y,Z, fft, back, FFTW_BACKWARD, FFTW_ESTIMATE);

  printf("plan built\n");

  // #pragma omp parallel for
  for (z=0;z<Z;z++) {
	for (y=0;y<Y;y++) {
	  for (x=0;x<X;x++) {
		index = z*Y*X + y*X + x;
		  orig[index] = 1.0*myrand() + 0.5;
		  if(index==0)
		  	printf("First is: %d\n", (long)orig[index]);
      sum += orig[index];
          //          printf("At %d we have %f %f (from %f)\n", index, orig[index][0], orig[index][1], 1.0*myrand() + 0.5);
	  }
	}
  }

  printf("values set\n");

  fftw_execute(p1);

  /*
  isys[0] = 1;
  zzfft3d(0, X, Y, Z, 1.0f, (complex *)sgdata, X, Y, (complex *)sgfftdata, X, Y,
		  table, work, isys);
  zzfft3d(1, X, Y, Z, 1.0f, (complex *)sgdata, X, Y, (complex *)sgfftdata, X, Y,
          table, work, isys);
  */
  printf("Forward FFT complete\n");


  fftw_execute(p2);
  //  csfft3d(-1, X, Y, Z, 1.0f/(1.0*X*Y*Z),(complex *) sgfftdata, X2, Y, 
  //		  (float *) sgdata, X+2, Y,table, work, isys);

   
  printf("Backward FFT complete!\n");

  //#pragma omp parallel for
  for (z=0;z<Z;z++) {
		for (y=0;y<Y;y++) {
	  	for (x=0;x<X;x++) {
				index = z*Y*X + y*X + x;
        back[index] /= X*Y*Z;
        // back[index] /= X*Y*Z;

				if( (fabs( (orig[index]-back[index])/(orig[index])) > precision) && orig[index]) {
					printf("Values are %d : %f %f\n", (int)index, orig[index],back[index]);
		  		total++;
				}
	  	}
		}
  }

  printf("[Double] Total number of voxels less than precision (%f %%): %d\n", precision*100, total);
  printf("[Double] Sum is: %ld, making a mean of %f\n", (long)(sum), (long)(sum)/((float)X*Y*Z));
/*
  for (z=0;z<Z;z++) {
    for (y=0;y<Y;y++) {
      for (x=0;x<X;x++) {
  		index = z*Y*X + y*X + x;
        //        printf("Values are %d : %f %f\n", (int)index, orig[index][0],back[index][0]);
  	  }
  	}
  }
  */
  //printf("End expected at %d\n", Z*Y*X2P1-1);

  fftw_destroy_plan(p1);
  fftw_destroy_plan(p2);
  fftw_free(orig); fftw_free(fft); fftw_free(back);
  // fftw_cleanup_threads();


  return 1;


}

int singletest() {

  //float sgdata[Z][Y][X];
  //complex sgfftdata[Z][Y][X/2];
  //float table[(Z + 256) + (2*Y + 256) + (2*X + 256)];
  //float work[Z + 4*X];
  fdrf_complex* orig = NULL;
  fdrf_complex* fft = NULL;
  fdrf_complex* back = NULL;
  fftwf_plan p1,p2;
  long x,y,z;
  long index;
  long total = 0;

  orig = (fdrf_complex *) fftw_malloc( X * Y * Z * sizeof(fdrf_complex));
  fft  = (fdrf_complex *) fftw_malloc( X * Y * Z * sizeof(fdrf_complex));
  back = (fdrf_complex *) fftw_malloc( X * Y * Z * sizeof(fdrf_complex));

  /* Build the plan */
  p1 = fftwf_plan_dft_3d(X,Y,Z, orig, fft, FFTW_FORWARD, FFTW_ESTIMATE);
  p2 = fftwf_plan_dft_3d(X,Y,Z, fft, back, FFTW_BACKWARD, FFTW_ESTIMATE);

  printf("plan built\n");

  // #pragma omp parallel for
  for (z=0;z<Z;z++) {
	for (y=0;y<Y;y++) {
	  for (x=0;x<X;x++) {
		index = z*Y*X + y*X + x;
		  orig[index] = 1.0*myrand() + 0.5;
          //          printf("At %d we have %f %f (from %f)\n", index, orig[index][0], orig[index][1], 1.0*myrand() + 0.5);
	  }
	}
  }

  printf("values set\n");

  fftwf_execute(p1);

  /*
  isys[0] = 1;
  zzfft3d(0, X, Y, Z, 1.0f, (complex *)sgdata, X, Y, (complex *)sgfftdata, X, Y,
		  table, work, isys);
  zzfft3d(1, X, Y, Z, 1.0f, (complex *)sgdata, X, Y, (complex *)sgfftdata, X, Y,
          table, work, isys);
  */
  printf("Forward FFT complete\n");


  fftwf_execute(p2);
  //  csfft3d(-1, X, Y, Z, 1.0f/(1.0*X*Y*Z),(complex *) sgfftdata, X2, Y, 
  //		  (float *) sgdata, X+2, Y,table, work, isys);

   
  printf("Backward FFT complete!\n");

  //#pragma omp parallel for
  for (z=0;z<Z;z++) {
		for (y=0;y<Y;y++) {
	  	for (x=0;x<X;x++) {
				index = z*Y*X + y*X + x;
				back[index] /= X*Y*Z;
				if( (fabs( (orig[index]-back[index])/(orig[index])) > precision) && orig[index]) {
			  	// printf("Values are %d : %f %f\n", (int)index, orig[index][0],back[index][0]);
			  	total++;
				}
	 		}
		}
  }

  printf("[Single] Total number of voxels less than precision (%f %%): %d\n", precision*100, total);

/*
  for (z=0;z<Z;z++) {
    for (y=0;y<Y;y++) {
      for (x=0;x<X;x++) {
  		index = z*Y*X + y*X + x;
        //        printf("Values are %d : %f %f\n", (int)index, orig[index][0],back[index][0]);
  	  }
  	}
  }
  */
  //printf("End expected at %d\n", Z*Y*X2P1-1);

  fftwf_destroy_plan(p1);
  fftwf_destroy_plan(p2);
  fftw_free(orig); fftw_free(fft); fftw_free(back);


  return 1;


}


int main() {

	// fftw_init_threads();
  // printf("**************** DOUBLE PRECISION TEST ********************\n");
  // doubletest();
  printf("*********************** 2D TEST ***************************\n");
  twodtest();
  // printf("**************** SINGLE PRECISION TEST ********************\n");
  // singletest();


}
