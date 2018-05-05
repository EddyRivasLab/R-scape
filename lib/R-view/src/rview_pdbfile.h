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
  int64_t  resnum;
  int64_t  dsq;
  
  ATOM    *atom;
  
} RES;

/* A chain has a sequence (SEQRES in a .pdb file)
 *
 * It also includes a collection of atoms per residue
 * The atom description includes a residue number (resnum), but
 * that number does not have to correspond to its position
 * in the SEQRES.
 *
 * This is a  huge pain. 
 *
 * We need to save both when reading a pdb file, and then
 * provide a map between the two:
 *
 *  atom_map[0..nr] in [0..L]
 *
 * Most methods don't bother producing the actual SEQRES or
 * giving the correspondence. They are happy by just reporting
 * the resnum in the atom description. 
 * That makes the results unmapable to any sequence
 *
 */
typedef struct chain_s {
  int64_t   L;
  ESL_DSQ  *seq;
  
  int64_t   nr;
  RES      *res;
  
  int      *atom_map;
} CHAIN;


#endif /*RVIEW_PDBFILE_INCLUDED*/
 
