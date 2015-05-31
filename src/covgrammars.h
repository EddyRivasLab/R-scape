/* covgrammarsr.h
 *
 *   
*/
#ifndef COVGRAMMARS_INCLUDED
#define COVGRAMMARS_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"


typedef float SCVAL;

typedef struct {
  int       L;    /* msa length */
  SCVAL  **dp;   /* L * L triangular DP matrix */
} GMX;

typedef struct {
} BGR;

typedef struct {
} G6;

typedef struct {
} BGR_MX;

typedef struct {
  GRM *S;
  GRM *L;
  GRM *F;
} G6_MX;

extern GMX  *GMX_Create(int L);
extern void  GMX_Destroy(GMX *gmx);
extern void  GMX_Dump(FILE *fp, GMX *gmx);


#endif
