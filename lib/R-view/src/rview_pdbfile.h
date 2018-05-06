/* rview_pddfile - functions to read a pdbx/mmcif file
 *
 */
#ifndef RVIEW_PDBFILE_INCLUDED
#define RVIEW_PDBFILE_INCLUDED

#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_sq.h"

typedef enum {
  ATYPE_ATOM   = 0,
  ATYPE_HETATM = 1,
  ATYPE_HOH    = 2,
  ATYPE_WAT    = 3,
  ATYPE_NONE   = 4,
} ATOMTYPE;

/* ATOM structure:
 *
 *  in pdbx/mmcif there are 21 entries per ATOM
 *
 *   ATOM  13786 C CA  . THR H  4 213 ? -15.552  -49.369 93.150  1.00 100.00 ? 209 THR L CA  1
 *
 * 
 *  0  _atom_site.group_PDB 
 *  1  _atom_site.id 
 *  2  _atom_site.type_symbol 
 *  3  _atom_site.label_atom_id 
 *  4  _atom_site.label_alt_id 
 *  5  _atom_site.label_comp_id 
 *  6  _atom_site.label_asym_id 
 *  7  _atom_site.label_entity_id 
 *  8  _atom_site.label_seq_id 
 *  9  _atom_site.pdbx_PDB_ins_code 
 *  10 _atom_site.Cartn_x 
 *  11 _atom_site.Cartn_y 
 *  12 _atom_site.Cartn_z 
 *  13 _atom_site.occupancy 
 *  14 _atom_site.B_iso_or_equiv 
 *  15 _atom_site.pdbx_formal_charge 
 *  16 _atom_site.auth_seq_id 
 *  17 _atom_site.auth_comp_id 
 *  18 _atom_site.auth_asym_id 
 *  19 _atom_site.auth_atom_id 
 *  20 _atom_site.pdbx_PDB_model_num 
 *
 */
typedef struct atom_s {
  char     *chain;   // chain
  ATOMTYPE  type;    // atom type ATOMT/HETATOMT
  ESL_DSQ   reschr;  // the residue 
  
  int64_t   idx;     // atom idx in the atom description

  /* the cartesian coordenates of the atom */
  double    x;       
  double    y;
  double    z;

} ATOM;

/* A residue: RES
 *   
 *      RES includes:
 *
 *           char *chain: the chain it belongs to
 */
typedef struct residue_s {
  char     *chain;         // chain the residue belongs to

  ESL_DSQ   reschr;        // the digitized character of RES
  int64_t   resnum;        // the number of RES in the atom description
 
  int64_t   from_atomidx;  // idx of the first ATOM associated to RES
  int64_t   na;            // number of atoms in RES
  ATOM     *atom;          // a ATOM structure for each atom
  
} RES;

/* A chain: CHAIN
 *
 *          CHAIN includes:
 ]*
 *                ESL_DSQ  *seq: the chain sequence. It can be found at:   
 *                                                      SEQRES in a .pdb file 
 *                                                      _entity_poly.pdbx_seq_one_letter_code_can in a .cif file
 *
 *                RES      *res: a collection of atoms per residue
 *
 * The atom description includes a residue number (resnum), but
 * that resnum does necessarily to correspond to its position inSEQRES.
 *
 * This is a HUGE pain. 
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
struct chain_s {
  int64_t   L;         // length of the seq
  ESL_SQ   *seq;       // The seq are read in the pdbfile (SEQRES)
  
  int64_t   nr;        // number of residues in the atom descripton
  RES      *res;       // a RES structure for each residue
  ESL_SQ   *resseq;    // The nr-long sequence with all residues
 
  int      *atom_map;  // atom_map[0..nr-1] with values in [0..L-1]
};

extern int  rview_ReadPDBxfile(char *pdbxfile, int *ret_nchain, struct chain_s **ret_chain, char *errbuf, int verbose);
extern void rview_atom_Init(ATOM *atom);
extern void rview_atom_Destroy(ATOM *atom);
extern int  rview_res_AddAtom(RES *res);
extern void rview_res_Init(RES *ret_res);
extern void rview_res_Destroy(RES *res);
extern int  rview_chain_AddRes(struct chain_s *chain);
extern void rview_chain_Init(struct chain_s *ret_chain);
extern void rview_chain_Destroy(struct chain_s *chain);

#endif /*RVIEW_PDBFILE_INCLUDED*/
 
