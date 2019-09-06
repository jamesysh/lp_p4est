#ifndef __TOOL_FN_H__
#define __TOOL_FN_H__
#include "initializer.h"

static  int
refine_init (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
    if(quadrant->level>4)

        return 0;
    else 
        return 1;
}


static  void        *

sc_array_index_begin (sc_array_t * arr)
{
  P4EST_ASSERT (arr != NULL);

  if (arr->elem_count == 0) {
    return NULL;
  }

  P4EST_ASSERT (arr->array != NULL);
  return (void *) arr->array;
}





#endif
