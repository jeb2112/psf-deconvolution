#include <fftw.h>
#include "math.h"
#include <stdio.h>
#include <stdlib.h>

fdr_complex *sgdata = NULL;
fdr_complex *sgfftdata = NULL;
float *table = NULL;
float *work = NULL;
int isys[2];

#define N 1024  

float gaussian(float x, float amp, float center, float lw) {

  return amp*exp(- (x-center)*(x-center) / (lw*lw) );

}

int fftshift(int n,fdr_complex *data) {

  fdr_complex *old_loc;
  fdr_complex *new_loc;
  fdr_complex temp;
  long k;

  for (k=0;k<n/2;k++) {
	old_loc = data + k;
	new_loc = data + k+n/2;
	temp = (*new_loc);
	(*new_loc) = (*old_loc);
	(*old_loc) = temp;
  }

  return 1;
}

int main() {
  fdr_complex x[N];
  fdr_complex y[N];
  float table[2*N+256] ;
  float work[2*N];
  int i;
  int isys[2];
  float sum = 0;

  for(i=0;i<N;i++) {
    x[i].re = gaussian(i*1.0f, 50.,N/2.0,N/10.0);
	x[i].im = 0;
  }

  for(i=0;i<N;i++) {
    fprintf(stdout,"%d = %f\n", i, x[i].re);
  }

  fftshift(N, x);

  for(i=0;i<N;i++) {
	fprintf(stdout,"%d = %f\n", i, x[i].re);
	sum += x[i].re;
  }
  
  fprintf(stdout,"Sum of %f\n", sum);
  
  isys[0] = 1;
  ccfft(0, N,1.0f,x,y,table,work,isys);
  ccfft(1, N,1.0f,x,y,table,work,isys);

  for(i=0;i<N;i++) {
    fprintf(stdout,"%d (re=%05.5f) (im=%05.5f) (abs=%05.5f)\n",i,y[i].re,y[i].im, sqrt(y[i].re*y[i].re+y[i].im*y[i].im));
  }

  return 1;

  

}
