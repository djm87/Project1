
#include <stdio.h>
#include <math.h>
#include "Project1.h"

// ----------------------------------------------------------------------
// main
//
// Solves for KdV Solution

int
main(int argc, char **argv)
{
  const int n = 100; // number of points
  const int tmax = 0.1, nprint = 25, a = 9, b = 16, c = 2;  
  double beg,end,dt,dx,dk;
  //int nmax; 
  
  dt = 0.1/pow(n,2); 
  dx = 20/n*2-20/n;
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
   VEC(x, i) = -10+20/n*i;
   
   if(i <= n/2) 
     VEC(k, i) = dk*i;
   else
     VEC(k, i) = -(n/2-1)*dk+dk*i;
   
   VEC(u,i) = 3*pow(a,2)*pow((1/cosh(0.5*a*VEC(x,i))),2);
  }

  // calculate scalar product and make sure it's correct
  //assert(x->n==y->n);

  beg =WTime();
  //#pragma omp parallel for
  //for (int i = 0; i < 10000; i++) {
  //   Interp_mid_add(x,y,interval, N);
  //}
  end = WTime();

  printf("Number of Points = %d\n", n);
  //for (int i = 0; i < x-> n; i++) {
    //printf("x(%2d) = %6.3f\t", i+1, VEC(x, i));
    //printf("y(%2d) = %6.3f\n", i+1, VEC(y, i));
   // fprintf(f, "%g %g\n", VEC(x, i), VEC(y, i));
  //}
  fclose(f);
  printf("Time = %.6f seconds\n", end-beg);
  //printf("Use './gnuplot.bash' to plot to file, then use 'eog Example1_plot.png'\n");

  // clean up
  vector_destroy(x);
  vector_destroy(k);
  vector_destroy(u);
   
  return 0;
}
