
#include <stdio.h>
#include <math.h>
#include <fftw.h>
#include <rfftw.h>
#include <complex.h>
#include "Project1.h"

#define PI (3.141592653589793)
#define tmax (1)
#define nprint (25)   
// ----------------------------------------------------------------------
// main
//
// Solves for KdV Solution
int
main(int argc, char **argv)
{
  const int n = 256; // number of points
  const int  a = 9,b = 16;// c = 2;  
  int imax,iplt,sprint; 
  double beg,beg2,end,dt,dx,dk,t,total_t,plotu[n*(nprint)][3],S_initializate,S_RK4_largest,S_Write_to_file ; 
  fftw_complex in[n], U[n], out[n],out2[n],
               ik3[n],g[n],E[n],E2[n],aa[n],bb[n],cc[n],dd[n];
  double t_cpu[10] = {0};
  fftw_plan plan_forward; 
  fftw_plan plan_backward;

  beg = WTime();

  dt = 0.1/(double)pow(n,2);
  //dt = 0.00001; 
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

    VEC(u,i) = 3*pow(a,2)*pow((1/cosh(0.5*a*VEC(x,i))),2)+
               3*pow(b,2)*pow((1/cosh(0.5*b*(VEC(x,i)+2))),2); //checked and good... 
 
    in[i].re = VEC(u,i);
    in[i].im = 0.0;

    ik3[i].im = pow(VEC(k,i),3); // checked and good...
  

  }

 fftw_one(plan_forward, in, U);

  end = WTime();
  t_cpu[0] = end-beg; 

 
  // Solver RK4
 
//printf("n=%d\n",nmax);

for (int i = 1; i<=imax ; i++) { //nmax; i++) { 
    beg2 = WTime();
    t=i*dt;
    
     for (int j = 0; j < x->n; j++)  { 
        g[j].im=-0.5*dt*VEC(k,j); //checked and good...
       E[j].re= cos(dt*ik3[j].im/2);
       E[j].im= sin(dt*ik3[j].im/2);

       E2[j].re = E[j].re*E[j].re-E[j].im*E[j].im;
       E2[j].im = 2*E[j].re*E[j].im; //checked and good...
     } 
      end = WTime();
  t_cpu[1] += end-beg2; 
     // a
     beg2 = WTime();
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
           end = WTime();
  t_cpu[2] += end-beg2; 
     // b
      beg2 = WTime();
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
                 end = WTime();
  t_cpu[3] += end-beg2; 
      // c
      beg2 = WTime();
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
               end = WTime();
  t_cpu[4] += end-beg2; 
      // d
      beg2 = WTime(); 
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
                 end = WTime();
  t_cpu[5] += end-beg2; 
      // Compute U 
      beg2 = WTime();
      for (int j = 0; j < x->n; j++) {
      U[j].re = E2[j].re*U[j].re - E2[j].im*U[j].im + 
               (E2[j].re*aa[j].re - E2[j].im*aa[j].im +
                2*(E[j].re*(bb[j].re+cc[j].re) - E[j].im*(bb[j].im+cc[j].im))  +
                dd[j].re)/6; 
        
      U[j].im = E2[j].re*U[j].im + E2[j].im*U[j].re + 
               (E2[j].re*aa[j].im + E2[j].im*aa[j].re + 
                2*(E[j].re*(bb[j].im+cc[j].im) + E2[j].im*(bb[j].re+cc[j].re)) +
                dd[j].im)/6;
      in[j].re = U[j].re;
      in[j].im = U[j].im;
      }  
           end = WTime();
  t_cpu[6] += end-beg2; 
      // Store data to print to file
      if(fmod(i,iplt) == 0){
      beg2 = WTime();
        fftw_one(plan_backward, in, out);
	//          for (int j = 0; j < x->n; j++) {
   // printf("i = %3d aa= %12.5f %12.5f \n",j, U[j].re/n ,U[j].im/n);
    //  } 
    //   printf("t = %f\n",t); 
 //getchar();
         for (int j = 0; j < x->n; j++) {
                 plotu[j+sprint][0] = t;
                  plotu[j+sprint][1] = VEC(x,j);
		  plotu[j+sprint][2] = out[j].re/n;
         }
      sprint+=n;
                end = WTime();
  t_cpu[7] += end-beg2; 
     }

   
  }

  //	          for (int j = 0; j < x->n; j++) {
   // printf("i = %3d aa= %12.5f %12.5f \n",j, U[j].re ,U[j].im);
    //  } 
    //   printf("t = %f\n",t); 
 //getchar();
  beg = WTime();
  FILE *f = fopen("Test1024.asc","w");
  for (int i = 0; i < n*(nprint); i++) {
    fprintf(f, "%0.3f %0.3f %0.3f\n", plotu[i][0],plotu[i][1],plotu[i][2]);
 }
  fclose(f);
  end = WTime();
  t_cpu[8] = end-beg;
  total_t = 0.0;
 for (int i = 0; i < 8; i++) {
  total_t+=t_cpu[i];
  }
 t_cpu[9] = total_t-t_cpu[8]-t_cpu[0];

  // See what effect parallization will have in RK4 loop,Initialize, and Write to file assuming 8 cores
  S_initializate = 1/(1-t_cpu[0]/total_t+ t_cpu[0]/total_t/8);
  S_RK4_largest = 1/(1-t_cpu[1]/imax/total_t+ t_cpu[1]/imax/total_t/8);
  S_Write_to_file = 1/(1-t_cpu[8]/total_t+ t_cpu[8]/total_t/8);

  printf(" Profiling (N = %d, dt = %g, T = %0.3f)\n",n,dt,tmax);
  printf("===================================================================\n");
  printf("Initialize            = %.6f sec                 %%Total = %5.3f\n", t_cpu[0],t_cpu[0]/total_t*100);
  printf("RK4 Preliminary       = %.6f sec  %%RK4 = %0.4f  %%Total = %5.3f\n", t_cpu[1]/imax,t_cpu[1]/imax/t_cpu[8]*100,t_cpu[1]/imax/total_t*100);
  printf("RK4 a (per loop)      = %.6f sec  %%RK4 = %0.4f  %%Total = %5.3f\n", t_cpu[2]/imax,t_cpu[2]/imax/t_cpu[8]*100,t_cpu[2]/imax/total_t*100);
  printf("RK4 b (per loop)      = %.6f sec  %%RK4 = %0.4f  %%Total = %5.3f\n", t_cpu[3]/imax,t_cpu[3]/imax/t_cpu[8]*100,t_cpu[3]/imax/total_t*100);
  printf("RK4 c (per loop)      = %.6f sec  %%RK4 = %0.4f  %%Total = %5.3f\n", t_cpu[4]/imax,t_cpu[4]/imax/t_cpu[8]*100,t_cpu[4]/imax/total_t*100);
  printf("RK4 d (per loop)      = %.6f sec  %%RK4 = %0.4f  %%Total = %5.3f\n", t_cpu[5]/imax,t_cpu[5]/imax/t_cpu[8]*100,t_cpu[5]/imax/total_t*100);
  printf("RK4 U (per loop)      = %.6f sec  %%RK4 = %0.4f  %%Total = %5.3f\n", t_cpu[6]/imax,t_cpu[6]/imax/t_cpu[8]*100,t_cpu[6]/imax/total_t*100);
  printf("RK4 store (per print) = %.6f sec  %%RK4 = %0.4f  %%Total = %5.3f\n", t_cpu[7]/nprint,t_cpu[7]/nprint/t_cpu[8]*100,t_cpu[7]/imax/total_t*100);   
  printf("Write to file         = %.6f sec                 %%Total = %5.3f\n", t_cpu[8],t_cpu[8]/total_t*100);
  printf("===================================================================\n");
  printf("RK4 total             = %.6f sec  %%Total = %0.3f\n", t_cpu[9],t_cpu[9]/total_t*100);
  printf("Wall time             = %.6f sec                 \n\n", total_t);

  printf(" Max Potential Speedup for parallel implimentation of 8 cores\n");
  printf("=================================================================\n");
  printf("Initialization        = %0.4f                  \n", S_initializate);
  printf("Largest part of RK4   = %0.4f                  \n", S_RK4_largest);
  printf("Write to File         = %0.4f                  \n", S_Write_to_file);
  // clean up
  vector_destroy(x);
  vector_destroy(k);
  vector_destroy(u);
  fftw_destroy_plan( plan_forward );
  fftw_destroy_plan( plan_backward );

  return 0;
}
