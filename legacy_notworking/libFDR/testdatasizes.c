

#include "fdr2.h"

int main()
{
  printf("int : %d\n", sizeof(int));
  printf("short : %d\n", sizeof(short));
  printf("long : %d\n", sizeof(long));
  printf("float : %d\n", sizeof(float));
  printf("double : %d\n", sizeof(double));
  printf("fftw_complex : %d\n", sizeof(fftw_complex));
  printf("fftwf_complex : %d\n", sizeof(fftwf_complex));
  printf("fdr_complex : %d\n", sizeof(fdr_complex));
  printf("mifcomplex : %d\n", sizeof(mifcomplex_t));
  printf("midcomplex : %d\n", sizeof(midcomplex_t));
  
// 	  cout << "Integer types:\n";
// 			    printTypeSize<char>("char");
// 	      printTypeSize<short>("short");
// 	        printTypeSize<int>("int");
// 		  printTypeSize<long>("long");
// 		    printTypeSize<long long>("long long");
// 
// 		      cout << "\nPointers:\n";
// 		        printTypeSize<void*>("void*");
// 
// 			  cout << "\nFloating point types:\n";
// 			    printTypeSize<float>("float");
// 			      printTypeSize<double>("double");
// 			        printTypeSize<long double>("long double");
// 
// 				  cout << "\nSizes from stddef.h:\n";
// 				    printTypeSize<size_t>("size_t");
// 				      printTypeSize<ptrdiff_t>("ptrdiff_t");

				        return 0;
}

