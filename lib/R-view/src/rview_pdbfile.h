/* rview_pddfile - functions to read a pdbx/mmcif file
 *
 */
#ifndef RVIEW_PDBFILE_INCLUDED
#define RVIEW_PDBFILE_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"

#include "rview_pdbfile.h"


typedef enum {
  ATOM       = 0,
  HETATOM    = 1,
} RESTYPE;

typedef struct atom_s {
  char     *chain;
  RESTYPE   type;
  ESL_DSQ   resdsq;
  int64_t   idx;

  double    x;
  double    y;
  double    z;
} ATOM;

typedef struct residue_s {
  char    *chain;
  int64_t  na;
  int64_t  from_atom;
  int64_t  num;
  int64_t  dsq;
  
  ATOM    *atom;
 
} RES;

typedef struct chain_s {
  int64_t  nr;
  RES     *res;
} CHAIN;


#endif /*RVIEW_PDBFILE_INCLUDED*/
 
