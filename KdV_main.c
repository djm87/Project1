
#include <stdio.h>
#include <math.h>
#include <fftw.h>
#include "Project1.h"

// ----------------------------------------------------------------------
// main
//
// Solves for KdV Solution

int
main(int argc, char **argv)
{
  const int n = 512; // number of points
  const int tmax = 0.1, nprint = 25, a = 9, b = 16, c = 2;  
  double beg,end,dt,dx,dk,t_init,t_rk4;
  //fftw_complex in[n], out[n]
  //fftw_plan p; 
  
   
  beg = WTime();

  dt = 0.1/pow(n,2); 
  dx = 20/n;
  dk = 2*pi/(n*dx);
  //nmax = round(tmax/dt); 
  //nplt = floor(tmax/nprint/dt);
  
  // init vectors
  struct vector *x = vector_create(n);
  struct vector *k = vector_create(n);
  struct vector *u = vector_create(n);

  FILE *f = fopen("Example1.asc","w");

  // initialize values 
 for (int i = 0; i < x->n; i++) {
   VEC(x, i) = -10+20.0/(n-1)*i; //checked and good... 
   
   if(i <= n/2) 
     VEC(k, i) = dk*i;
   else
     VEC(k, i) = -(n/2-1)*dk+dk*i;
   
   //VEC(u,i) = 3*pow(a,2)*pow((1/cosh(0.5*a*VEC(x,i))),2); //checked and good... 
   //printf("u(%2d) = %6.4f\t", i+1, VEC(u, i));
  }
 end = WTime();
  t_init = end-beg; 
  printf("Number of Points = %d\n", n);
  //p = fftw_create_plan(n,FFTW_FORWARD,FFTW_ESTIMATE);
  //fftw_one(p,in,out); 
  //assert(x->n==y->n);

  // Solver RK4
  beg = WTime();
//
  end = WTime();
  t_rk4 = end-beg;

  fclose(f);
  printf("Time Initialize = %.6f seconds\n", t_init);
  printf("Time RK4 = %.6f seconds\n", t_rk4);

  // clean up
  vector_destroy(x);
  vector_destroy(k);
  vector_destroy(u);
  //fftw_destroy_plan(p) 
  return 0;
}
