
#include <stdio.h>
#include <math.h>
#include <fftw.h>
#include <rfftw.h>
#include <complex.h>
#include "Project1.h"

#define PI (3.141592653589793)
#define tmax (0.1)
#define nprint (25)   
// ----------------------------------------------------------------------
// main
//
// Solves for KdV Solution

int
main(int argc, char **argv)
{
  const int n = 512; // number of points
  const int  a = 9;// b = 16, c = 2;  
  int imax,iplt,sprint; 
  double beg,end,dt,dx,dk,t,t_cpu[3],total_t, plotu[n*nprint][3] ; 
  fftw_complex in[n], U[n], out[n],out2[n],
               ik3[n],g[n],E[n],E2[n],aa[n],bb[n],cc[n],dd[n];

  fftw_plan plan_forward; 
  fftw_plan plan_backward;

  beg = WTime();

  dt = 0.1/(double)pow(n,2); 
  dx = 20/((double)n-1);
  dk = 2*PI/(dx*(double)n);
  imax = round(tmax/dt);
  iplt = floor(tmax/nprint/dt);
  sprint=0;
  // init vectors
  struct vector *x = vector_create(n);
  struct vector *k = vector_create(n);
  struct vector *u = vector_create(n);
  
  
  plan_forward = fftw_create_plan(n,FFTW_FORWARD,FFTW_ESTIMATE);
  plan_backward = fftw_create_plan(n,FFTW_BACKWARD,FFTW_ESTIMATE);


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

    ik3[i].im = pow(VEC(k,i),3); // checked and good...
  

  }

 fftw_one(plan_forward, in, U);

  end = WTime();
  t_cpu[0] = end-beg; 

 
  // Solver RK4
  beg = WTime();
//printf("n=%d\n",nmax);

for (int i = 0; i<imax ; i++) { //nmax; i++) { 
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

      // c
      for (int j = 0; j < x->n; j++) {
      in[j].re = E[j].re*U[j].re - E[j].im*U[j].im+bb[j].re/2 ;
      in[j].im = E[j].re*U[j].im + E[j].im*U[j].re+bb[j].im/2 ; // checked and good...
      }   
  
      fftw_one(plan_backward, in, out2);
     
      for (int j = 0; j < x->n; j++) {
      in[j].re = (out2[j].re/n)*(out2[j].re/n); 
      in[j].im = 0.0; // checked and good...
      }

      fftw_one(plan_forward, in, out);
      for (int j = 0; j < x->n; j++) {
      cc[j].re = -g[j].im*out[j].im;
      cc[j].im = g[j].im*out[j].re;
      }  
    
      // d 
      for (int j = 0; j < x->n; j++) {
      in[j].re = E2[j].re*U[j].re - E2[j].im*U[j].im + E[j].re*cc[j].re - E[j].im*cc[j].im;
      in[j].im = E2[j].re*U[j].im + E2[j].im*U[j].re + E[j].re*cc[j].im + E[j].im*cc[j].re ; // checked and good...
      }   
  
      fftw_one(plan_backward, in, out2);
     
      for (int j = 0; j < x->n; j++) {
      in[j].re = (out2[j].re/n)*(out2[j].re/n); 
      in[j].im = 0.0; // checked and good...
      }

      fftw_one(plan_forward, in, out);
      for (int j = 0; j < x->n; j++) {
      dd[j].re = -g[j].im*out[j].im;
      dd[j].im = g[j].im*out[j].re;
      }  
      
      // Compute U 
      for (int j = 0; j < x->n; j++) {
      U[j].re = E2[j].re*U[j].re - E2[j].im*U[j].im + 
               (E2[j].re*aa[j].re - E2[j].im*aa[j].im +
                2*(E[j].re*(bb[j].re+cc[j].re) - E[j].im*(bb[j].im+cc[j].im))  +
                dd[j].re)/6; 
        
      U[j].im = E2[j].re*U[j].im + E2[j].im*U[j].re + 
               (E2[j].re*aa[j].im + E2[j].im*aa[j].re + 
                2*(E[j].re*(bb[j].im+cc[j].im) + E2[j].im*(bb[j].re+cc[j].re)) +
                dd[j].im)/6;
      }  
      // Store data to print to file
      if(fmod(i,iplt) == 0){
         fftw_one(plan_backward, U, out);
	     
         for (int j = 0; j < x->n; j++) {
                  plotu[j+sprint][0] = t;
                  plotu[j+sprint][1] = VEC(x,j);
		  plotu[j+sprint][2] = out[j].re;
         }
      sprint+=n;
      }
  }
 //for (int i = 0; i < x->n; i++) {
 //   printf("i = %3d aa= %12.4g %12.4g \n",i, U[i].re ,U[i].im);
// }
  end = WTime();
  t_cpu[1] = end-beg;

  
  beg = WTime();
  FILE *f = fopen("Test1.asc","w");
  for (int i = 0; i < n*nprint; i++) {
    fprintf(f, "%0.3f %0.3f %0.3f\n", plotu[i][0],plotu[i][1],plotu[i][2]);
  }
  fclose(f);
  end = WTime();
  t_cpu[2] = end-beg;

  total_t = t_cpu[0]+t_cpu[1]+t_cpu[2];

  printf("                      Profiling\n");
  printf("=====================================================\n");
  printf("Time Write to file  = %.6f sec  %%Total = %0.3f\n", t_cpu[2],t_cpu[2]/total_t*100);
  printf("Time Initialize     = %.6f sec  %%Total = %0.3f\n", t_cpu[0],t_cpu[0]/total_t*100);
  printf("Time RK4            = %.6f sec  %%Total = %0.3f\n", t_cpu[1],t_cpu[1]/total_t*100);
  
  // clean up
  vector_destroy(x);
  vector_destroy(k);
  vector_destroy(u);
  fftw_destroy_plan( plan_forward );
  fftw_destroy_plan( plan_backward );

  return 0;
}
