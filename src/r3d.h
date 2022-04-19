/* r3d.h
 *
 *   
*/
#ifndef R3D_INCLUDED
#define R3D_INCLUDED

#include <stdio.h>		/* FILE */


#include "easel.h"
#include "esl_random.h"
#include "esl_sq.h"

#include "e2_profilesq.h"
#include "covgrammars.h"
#include "correlators.h"
#include "structure.h"

typedef struct {
} R3D_HL;

typedef struct {
} R3D_HLparam;

typedef struct {
} R3D_BL;

typedef struct {
} R3D_BLparam;

typedef struct {
} R3D_IL;
typedef struct {
} R3D_ILparam;

typedef struct {
  int     nHL;
  R3D_HL *HL;
  
  int     nBL;
  R3D_BL *BL;
  
  int     nIL;
  R3D_IL *IL;
} R3D;

typedef struct {
  R3D_HLparam *HLp;
  R3D_BLparam *BLp;
  R3D_ILparam *ILp;
} R3Dparam;

typedef struct {
} R3D_HLMX;
 
typedef struct {
} R3D_BLMX;

typedef struct {
} R3D_ILMX;

typedef struct {
  R3D_HLMX *HLmx;
  R3D_BLMX *BLmx;
  R3D_ILMX *ILmx;
} R3D_MX;

const R3D_HLparam R3D_HL_PRELOADS = {
};
const R3D_BLparam R3D_BL_PRELOADS = {
};
const R3D_ILparam R3D_IL_PRELOADS = {
};

extern int R3D_GetParam (R3Dparam  **ret_r3dp, char *errbuf, int verbose);

#endif
