
#include <stdio.h>
#include <math.h>
#include <fftw.h>
#include <rfftw.h>
#include <complex.h>
#include "Project1.h"

#define PI (3.141592653589793)
// ----------------------------------------------------------------------
// main
//
// Solves for KdV Solution

int
main(int argc, char **argv)
{
  const int n = 64; // number of points
  const int nprint = 25, a = 9;// b = 16, c = 2;  
  int nmax,nplt; 
  const double tmax = 0.01;
  double beg,end,dt,dx,dk,t_init,t_rk4,t; 
  fftw_complex in[n], U[n], in2[n], out[n],out2[n],ik3[n],g[n],E[n],E2[n],aa[n],bb[n],cc[n],dd[n];
  fftw_real  inr[n], outr[n], inr2[n], outr2[n];

  fftw_plan plan_forward; 
  fftw_plan plan_backward;
  fftw_plan  plan_r_to_c;
  fftw_plan  plan_c_to_r;

  beg = WTime();

  dt = 0.1/(double)pow(n,2); 
  dx = 20/((double)n-1);
  dk = 2*PI/(dx*(double)n);
  nmax = round(tmax/dt);
  nplt = floor(tmax/nprint/dt);
  
  // init vectors
  struct vector *x = vector_create(n);
  struct vector *k = vector_create(n);
  struct vector *u = vector_create(n);
  
  plan_forward = fftw_create_plan(n,FFTW_FORWARD,FFTW_ESTIMATE);
  plan_backward = fftw_create_plan(n,FFTW_BACKWARD,FFTW_ESTIMATE);
  
  plan_r_to_c = rfftw_create_plan(n,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);
  plan_c_to_r = rfftw_create_plan(n,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE);

  FILE *f = fopen("Example1.asc","w");

  // initialize values 
  for (int i = 0; i < x->n; i++) {
    VEC(x, i) = -10+20.0/(n-1)*i; //checked and good...   
    if(i <= n/2){
      VEC(k, i) = dk*i;}
    else
      VEC(k, i) = -(n/2-1)*dk+(i-1-n/2)*dk;

    VEC(u,i) = 3*pow(a,2)*pow((1/cosh(0.5*a*VEC(x,i))),2); //checked and good... 
 
    in[i].re = VEC(u,i);
    in[i].im = 0.0;
    inr[i] = VEC(u,i); 

    ik3[i].im = pow(VEC(k,i),3); // checked and good...
  

  }

 fftw_one(plan_forward, in, U);

  end = WTime();
  t_init = end-beg; 

 
  // Solver RK4
  beg = WTime();
//printf("n=%d\n",nmax);

for (int i = 0; i<1 ; i++) { //nmax; i++) { 
    t=i*dt;
    
     for (int j = 0; j < x->n; j++)  { 
        g[j].im=-0.5*dt*VEC(k,j); //checked and good...
       E[j].re= cos(dt*ik3[j].im/2);
       E[j].im= sin(dt*ik3[j].im/2);

       E2[j].re = E[j].re*E[j].re-E[j].im*E[j].im;
       E2[j].im = 2*E[j].re*E[j].im; //checked and good...
     } 
     // a
     fftw_one(plan_backward, U, out);
     
      for (int j = 0; j < x->n; j++) {
      in[j].re = (out[j].re/n)*(out[j].re/n);
      in[j].im = 0.0; // checked and good...
      }

     fftw_one(plan_forward, in, out);
     for (int j = 0; j < x->n; j++) {
     aa[j].re = -g[j].im*out[j].im;
     aa[j].im = g[j].im*out[j].re;
     }
     
     // b
      for (int j = 0; j < x->n; j++) {
      in[j].re = E[j].re*(U[j].re+aa[j].re/2) - E[j].im*(U[j].im+aa[j].im/2) ;
      in[j].im = E[j].re*(U[j].im+aa[j].im/2)  + E[j].im*(U[j].re+aa[j].re/2) ; // checked and good...
      }   

      fftw_one(plan_backward, in, out2);
     
      for (int j = 0; j < x->n; j++) {
      in[j].re = (out2[j].re/n)*(out2[j].re/n);
      in[j].im = 0.0; // checked and good...
      }

      fftw_one(plan_forward, in, out);
      for (int j = 0; j < x->n; j++) {
      bb[j].re = -g[j].im*out[j].im;
      bb[j].im = g[j].im*out[j].re;
      }
  }
 for (int i = 0; i < x->n; i++) {
printf("i = %3d aa= %12.4g %12.4g \n",i, bb[i].re ,bb[i].im);
 }
  end = WTime();
  t_rk4 = end-beg;

  fclose(f);
  printf("Time Initialize = %.6f seconds\n", t_init);
  printf("Time RK4 = %.6f seconds\n", t_rk4);

  // clean up
  vector_destroy(x);
  vector_destroy(k);
  vector_destroy(u);
  fftw_destroy_plan( plan_forward );
  fftw_destroy_plan( plan_backward );
  fftw_destroy_plan(plan_r_to_c);
  fftw_destroy_plan(plan_c_to_r);  

  return 0;
}

//for (int i = 0; i < x->n; i++) {
  //in[i].im=out[i].im;
 // in[i].im=out[i].re; 
//}

//for (int i = 0; i < x->n; i++) {
//printf("i = %3d out= %12g %12g\n",i,out[i].re,out[i].im);
// }
