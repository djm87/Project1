
#include "Project1.h"
#include <fftw.h>
#include <assert.h>
#include <stdlib.h>

// ----------------------------------------------------------------------
// vector_create
//
// allocates and initializes a vector, setting all elements to zero

fftw_complex  
fftw_multiple(fftw_complex *a, fftw_complex *b, int n)
{
  fftw_complex *c[n];
  // This is simple for now, should add cases for real/complex cases
  for (int i = 0; i < n; i++) {
  {  
     c[i].re = a[i].re*b[i].re - a[i].im*b[i].im; 
     c[i].im = a[i].re*b[i].im + a[i].im*b[i].re;
  }
  return c;
}

fftw_complex 
fftw_add(fftw_complex *a, fftw_complex *b, int n)
{
  fftw_complex c[n];
  // add fftw complex
  for (int i = 0; i < n; i++) {
  {  
     c[i].re = a[i].re + b[i].re; 
     c[i].im = a[i].im + b[i].im;
  }
  return c;
}
