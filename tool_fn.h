#ifndef __TOOL_FN_H__
#define __TOOL_FN_H__


static int
refine_init (p4est_t * p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t * quadrant)
{
    if(quadrant->level>4)

        return 0;
    else 
        return 1;
}







#endif
