
#include "Project1.h"

#include <assert.h>
#include <stdlib.h>

// ----------------------------------------------------------------------
// vector_create
//
// allocates and initializes a vector, setting all elements to zero

struct vector *
vector_create(int n)
{
  // allocate vector struct
  struct vector *v = calloc(1, sizeof(*v));

  // remember how many elements the vector contains
  v->n = n;
  // allocate space for n elements and zero them.
  v->vals = calloc(n, sizeof(v->vals[0]));

  return v;
}

// ----------------------------------------------------------------------
// vector_destroy
//
// this function is called when we're done using a vector

void
vector_destroy(struct vector *v)
{
  free(v->vals);
  // not strictly needed, but will make sure that we crash if we
  // access the vector after we called vector_destroy()
  v->vals = NULL;
  free(v);
}

// ----------------------------------------------------------------------
// vector_is_equal
//
// checks two vectors for equality

bool
vector_is_equal(struct vector *v1, struct vector *v2)
{
  // make sure the vectors have the same length
  assert(v1->n == v2->n);

  for (int i = 0; i < v1->n; i++) {
    if (VEC(v1, i) != VEC(v2, i))
      return false;
  }
  return true;
}

