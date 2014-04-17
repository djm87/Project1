
// include guards: an idiom which ensures that it's okay to include
// the header more than once -- the 2nd time around, LINEAR_ALGEBRA_H will
// be defined already and so nothing will be included anymore.

#ifndef Project1_H
#define Project1_H




#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#define BOUNDSCHECK
# define pi 3.141592653589793238463
// ----------------------------------------------------------------------
// struct vector

struct vector {
  double *vals;
  int n;
};

#ifdef BOUNDSCHECK
#define VEC(v, i) (*({				\
	assert((i) >= 0 && (i) < (v)->n);	\
	&((v)->vals[(i)]);			\
      })) 
#else
#define VEC(v, i) ((v)->vals[i])
#endif

struct vector *vector_create(int n);
void vector_destroy(struct vector *v);
bool vector_is_equal(struct vector *v1, struct vector *v2);
void Interp_mid_add(struct vector *x, struct vector *y,const double interval, int N);
double trapz(int N);
double fnc(double x);
// ----------------------------------------------------------------------
// other useful stuff

#include <sys/time.h>

static inline double
WTime(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec / 1e6;
}

// ----------------------------------------------------------------------

#include <stdio.h>

#define HERE fprintf(stderr, "HERE at %s:%d (%s)\n", __FILE__, __LINE__, __FUNCTION__) 

#endif
