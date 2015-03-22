#ifndef BOYS_PARAM_H
#define BOYS_PARAM_H

// include parameters from here as well
#include "boys_shortgrid.h"

// These don't only apply to the short-range grid
// so sopy the definitions
#define BOYS_GRID_MAXN       BOYS_SHORTGRID_MAXN
#define BOYS_GRID_SPACE      BOYS_SHORTGRID_SPACE
#define BOYS_GRID_LOOKUPFAC  BOYS_SHORTGRID_LOOKUPFAC
#define BOYS_GRID_LOOKUPFAC2 BOYS_SHORTGRID_LOOKUPFAC2

#define BOYS_INTERP_ORDER 7

#endif
