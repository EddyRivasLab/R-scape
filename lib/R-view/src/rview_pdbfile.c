/* rview_pddfile - functions to read a pdbx/mmcif file
 * Contents:
 *
 * ER, Thu May  3 19:21:13 EDT 2018 [Verril Farm] 
 * SVN $Id:$
 */

#include "rview_config.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <math.h>
#include <float.h>
	
#include "easel.h"
#include "esl_buffer.h"
#include "esl_mem.h"
#include "esl_sq.h"

#include "rview_pdbfile.h"

static int       atmlabel_create  (int *ret_nal, char ***ret_atmlabel,    int verbose);
static int       sasymlabel_create(int *ret_nsl, char ***ret_sasymlabel,  int verbose);
static int       parse_atom             (PDBX *pdbx, char *p, esl_pos_t n, char *atomid,    int natf, int *useatmf,                        int verbose);
static int       parse_atomfield        (            char *p, esl_pos_t n, int nal, char **atmlabel, int *ret_atmidx, int **ret_useatmf,   int verbose);
static int       parse_pdbname          (            char *p, esl_pos_t n, char *entryid,   char **ret_pdbname,                            int verbose);
static int       parse_reslow           (            char *p, esl_pos_t n, char *rlow,      double *ret_resolution,                        int verbose);
static int       parse_reshigh          (            char *p, esl_pos_t n, char *rhigh,     double *ret_resolution,                        int verbose);
static int       parse_seqid            (PDBX *pdbx, char *p, esl_pos_t n, char *seq_id,    int64_t *ret_seqid,     int *ret_multiseq,     int verbose);
static int       parse_seqtype          (            char *p, esl_pos_t n, char *seq_type,  char **ret_seqtype,                            int verbose);
static int       parse_sequence_short   (            char *p, esl_pos_t n, char *seq_can,   char **ret_seq,                                int verbose);
static int       parse_sequence_long    (            char *p, esl_pos_t n,                  char **ret_seq,         esl_pos_t *ret_seqlen, int verbose);
static int       parse_sequence_end     (            char *p, esl_pos_t n, char *scomma,    int *ret_readnow,       struct chain_s *chain, int verbose);
static int       parse_sequence_hetero  (PDBX *pdbx, char *p, esl_pos_t n, int verbose);
static int       parse_multiseq         (            char *p, esl_pos_t n,                  char **ret_chunk,       esl_pos_t *ret_len,    int verbose);
static int       parse_add_allseqs      (PDBX *pdbx, int *ret_readnow, char *chunk, esl_pos_t len, int verbose);
static int       parse_sasymfield       (            char *p, esl_pos_t n, int nsl, char **sasymlabel, int *ret_nsaf, int **ret_usesasymf, int verbose);
static int       parse_sasym            (PDBX *pdbx, char *p, esl_pos_t n,                  int nsaf, int *usesasymf,                      int verbose);
static int       chunk_is_newsequence   (PDBX *pdbx, char *line, struct chain_s *chain, int *ret_isshortchain, int verbose);
static void      chunk_is_sequence_short(char *line, char **ret_seq, char **ret_chainname, int verbose);
static int       chunk_is_sequence_start(char *line, char **ret_seq, SEQSTATUS *ret_seqstatus);
static int       chunk_is_sequence_cont (char *line, char **ret_seq, SEQSTATUS *ret_seqstatus);
static int       chunk_is_sequence_end  (char *line,                 SEQSTATUS *ret_seqstatus);
static int       chunk_is_chainname(char *line, char **ret_name);
static int       res_oneletter(char **ret_x, char *aa);
static void      writep(char *p, esl_pos_t n);


/* rview_ReadPDBxfile()
 *
 */
int
rview_ReadPDBxfile(char *pdbxfile, PDBX **ret_pdbx, char *errbuf, int verbose)
{
  ESL_BUFFER      *bf      = NULL;
  PDBX            *pdbx    = NULL;
  struct chain_s  *chain   = NULL;
  char            *chunk   = NULL;
  char            *p;
  char             first[]         = "data_";
  char             entryid[]       = "_entry.id";
  char             seq_id[]        = "_entity_poly.entity_id";
  char             seq_type[]      = "_entity_poly.type";
  char             seq_can[]       = "_entity_poly.pdbx_seq_one_letter_code_can";
  char             seq_chainname[] = "_entity_poly.pdbx_strand_id";
  char             seq_targetid[]  = "_entity_poly.pdbx_target_identifier";
  char             seq_hetero[]    = "_entity_poly_seq.hetero";
  char             sasym[]         = "_struct_asym.";
  char             atomsite[]      = "_atom_site.";
  char             rhigh[]         = "_refine.ls_d_res_high";
  char             rlow[]          = "_refine.ls_d_res_low";
  char             isrna[]         = "polyribonucleotide";
  char             ispep[]         = "polypeptide";
  char             scomma[]        = ";";
  char             newcommand[]    = "_";
  char             endloop[]       = "#";
  char             atomid[]        = "ATOM";
  char             hatomid[]       = "HETATM";
  esl_pos_t        n;
  esl_pos_t        len = 0;
  int              multiseq;
  int              readnow;
  int              checkhetero = FALSE;
  int              nsl;
  int              nal;
  int              s;
  int              l;
  char           **sasymlabel = NULL;
  char           **atmlabel   = NULL;
  int             *useatmf    = NULL;
  int             *usesasymf  = NULL;
  int              natf = 0;
  int              nsaf = 0;
  int              issasym  = FALSE;
  int              verbose2 = TRUE;
  int              status;
  
  pdbx  = rview_pdbx_Create();  if (pdbx  == NULL) esl_fatal("could not create PDBX  structure");
  chain = rview_chain_Create(); if (chain == NULL) esl_fatal("could not create chain structure");

  // the struct_asym fields
  sasymlabel_create(&nsl, &sasymlabel, verbose);
  
  // the atom fields
  atmlabel_create(&nal, &atmlabel, verbose);
  
  // Read the pdbfile
  //
  // For each chain,
  // we are going to extract the sequence (chain->seq) and the information about all the residues and their atoms (chain->res)
  //
  status = esl_buffer_Open(pdbxfile, NULL, &bf);
  if      (status == eslENOTFOUND) esl_fatal("open failed: %s",     bf ? bf->errmsg : "(no other diagnostics available)");
  else if (status == eslFAIL)      esl_fatal("gzip -dc failed: %s", bf ? bf->errmsg : "(no other diagnostics available)");
  else if (status != eslOK)        esl_fatal("open failed with error code %d", status);
  
  // first line of file:
  //
  // data_pdbname
  // #
  //
  if ((status = esl_buffer_Get(bf, &p, &n)) != eslOK) goto ERROR; /* includes normal EOF */
  if (!esl_memstrpfx(p, n, first)) ESL_XFAIL(eslEINVAL, bf->errmsg, "Expected FASTA record to start with %s", first);
  
  while (( status = esl_buffer_GetLine(bf, &p, &n)) == eslOK) 
    {      
      // THE PDBNAME
      // # 
      // _entry.id   4RUM 
      // #
      //
      if (esl_memstrpfx(p, n, entryid)) parse_pdbname(p, n, entryid, &pdbx->pdbname, verbose2);
      
      // THE RESOLUTION
      //
      // _refine_hist.d_res_high                       2.6426 
      // _refine_hist.d_res_low                        33.867
      //
      if (esl_memstrpfx(p, n, rlow))     parse_reslow  (p, n, rlow,      &pdbx->res_low,  verbose2);
      if (esl_memstrpfx(p, n, rhigh))    parse_reshigh (p, n, rhigh,     &pdbx->res_high, verbose2);
      
      // Determine if there are multiple chains by looking at field
      // _entity_poly.entity_id
      //
      // if there is anumber after,  
      //
      // _entity_poly.entity_id                      1
      //
      //  there is only one sequence. Otherwise there will be an empty field as the sequences are
      //  reporterd all afterwarads together
      //
      // _entity_poly.entity_id                      
      //
      //
      if (esl_memstrpfx(p, n, seq_id)) parse_seqid(pdbx, p, n, seq_id, &chain->seqid, &multiseq, verbose2);

      // THE SEQUENCES
      // For the sequence, there are 8 fields, and
      // 2 possible formats (single-sequence vs multi-sequence),
      // and within each format, two sequence formats (short-seq vs long-seq)
      //
      // fields are:
      //
      // _entity_poly.type 
      // _entity_poly.nstd_linkage 
      // _entity_poly.nstd_monomer 
      // _entity_poly.pdbx_seq_one_letter_code 
      // _entity_poly.pdbx_seq_one_letter_code_can 
      // _entity_poly.pdbx_strand_id 
      // _entity_poly.pdbx_target_identifier 
      //
      // If SINGLE SEQUENCE, the value of the field seems to follow the field
      // for each sequence as,
      //
      // _entity_poly.entity_id                      1 
      // _entity_poly.type                           polyribonucleotide 
      // _entity_poly.nstd_linkage                   no 
      // _entity_poly.nstd_monomer                   yes 
      // _entity_poly.pdbx_seq_one_letter_code       
      // ;(GTP)GGAACUGAGCAGGCAAUGACCAGAGCGGUCAUGCAGCCGGGCUGCGAAAGCGGCAACAGAUGAUUACACGCACAU
      // CUGUGGGACAGUUCCCAC
      // ;
      // _entity_poly.pdbx_seq_one_letter_code_can   
      // ;GGGAACUGAGCAGGCAAUGACCAGAGCGGUCAUGCAGCCGGGCUGCGAAAGCGGCAACAGAUGAUUACACGCACAUCUGU
      // GGGACAGUUCCCAC
      // ;
      // _entity_poly.pdbx_strand_id                 A 
      // _entity_poly.pdbx_target_identifier         ? 
      // #
      //
      // but watch out for small sequences which come as
      //
      // _entity_poly.entity_id                      1 
      // _entity_poly.type                           polyribonucleotide 
      // _entity_poly.nstd_linkage                   no 
      // _entity_poly.nstd_monomer                   no 
      // _entity_poly.pdbx_seq_one_letter_code       GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG 
      // _entity_poly.pdbx_seq_one_letter_code_can   GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG 
      // _entity_poly.pdbx_strand_id                 A 
      // _entity_poly.pdbx_target_identifier         ? 
      // # 
      //
      //
      // If MULTISEQ, the format seems to be to put the field descriptions first, and then all the 8 fields
      // for each sequence as:
      //
      // _entity_poly.type 
      // _entity_poly.nstd_linkage 
      // _entity_poly.nstd_monomer 
      // _entity_poly.pdbx_seq_one_letter_code 
      // _entity_poly.pdbx_seq_one_letter_code_can 
      // _entity_poly.pdbx_strand_id 
      // _entity_poly.pdbx_target_identifier 
      // 1  'polypeptide(L)'   no no MKTDFNQKIEQLKEFIEECRRVWLVLKKPTKDEYLAVAKVTALGISLLGIIGYIIHVPATYIKGILK 
      // 1MKTDFNQKIEQLKEFIEECRRVWLVLKKPTKDEYLAVAKVTALGISLLGIIGYIIHVPATYIKGILK A7    ? 
      // 12  'polypeptide(L)'   no no MSKREETGLATSAGLIRYMDETFSKIRVKPEHVIGVTVAFVIIEAILTYGRF 
      // 1MSKREETGLATSAGLIRYMDETFSKIRVKPEHVIGVTVAFVIIEAILTYGRF A8    ? 
      // 13  'polypeptide(L)'   no no MARNKPLAKKLRLAKALKQNRRVPVWVIVKTNRRVLTHPKRRYWRRTKLKE 
      // 1MARNKPLAKKLRLAKALKQNRRVPVWVIVKTNRRVLTHPKRRYWRRTKLKE Af    ? 
      // 14  'polypeptide(L)'   no no 
      // 1;MNTLKMQRRIAAEILKCGENRIWIDPERIDDVASAITREDIKRLIKEGVIKKKPIKGQSRYRAKIRHEQKKKGRHRGPGS
      // 1RKGKKTARMGKKELWIKTIRALRKELRKLKAQKKIDRKTYRMLYIRAKGGQFKNKHQLYLFLEEHGLLKK
      // 1;
      // 1;MNTLKMQRRIAAEILKCGENRIWIDPERIDDVASAITREDIKRLIKEGVIKKKPIKGQSRYRAKIRHEQKKKGRHRGPGS
      // 1RKGKKTARMGKKELWIKTIRALRKELRKLKAQKKIDRKTYRMLYIRAKGGQFKNKHQLYLFLEEHGLLKK
      // 1;
      // 1AQ    ? 
      //
      // sequence format:
      // _entity_poly.pdbx_seq_one_letter_code 
      // _entity_poly.pdbx_seq_one_letter_code_can
      //
      // 4 different formats are used depending on whether the sequence is small and fits in one line or not
      //
      // (1) sort sequence:
      //
      // 1  'polypeptide(L)'   no no MKTDFNQKIEQLKEFIEECRRVWLVLKKPTKDEYLAVAKVTALGISLLGIIGYIIHVPATYIKGILK 
      // MKTDFNQKIEQLKEFIEECRRVWLVLKKPTKDEYLAVAKVTALGISLLGIIGYIIHVPATYIKGILK A7    ? 
      //
      // (2) long sequence
      //
      // 14  'polypeptide(L)'   no no 
      // ;MNTLKMQRRIAAEILKCGENRIWIDPERIDDVASAITREDIKRLIKEGVIKKKPIKGQSRYRAKIRHEQKKKGRHRGPGS
      // 1RKGKKTARMGKKELWIKTIRALRKELRKLKAQKKIDRKTYRMLYIRAKGGQFKNKHQLYLFLEEHGLLKK
      // ;
      // ;MNTLKMQRRIAAEILKCGENRIWIDPERIDDVASAITREDIKRLIKEGVIKKKPIKGQSRYRAKIRHEQKKKGRHRGPGS
      // 1RKGKKTARMGKKELWIKTIRALRKELRKLKAQKKIDRKTYRMLYIRAKGGQFKNKHQLYLFLEEHGLLKK
      // ;
      // 1AQ    ?
      //
      // (3) super short sequence
      //
      // 37 'polypeptide(L)'   no no MGQKWKLYEIKDGKVIRKNKFCPRCGPGVFMADHGDRWACGKCGYTEWKK MGQKWKLYEIKDGKVIRKNKFCPRCGPGVFMADHGDRWACGKCGYTEWKK
      // BY    ?
      //
      // (4) super-super short sequence
      //
      // 53 'polypeptide(L)'   no no MRWKWIKKRIRRLKRQRKKERGLI MRWKWIKKRIRRLKRQRKKERGLI Ah    ? 
      //
      //
      if (multiseq == FALSE) { // only one sequence
	// determine the type of sequence 
	if (esl_memstrpfx(p, n, seq_type)) parse_seqtype(p, n, seq_type, &chain->seqtype, verbose2);
	// determine when to read the canonical sequence 
	if (esl_memstrpfx(p, n, seq_can))  readnow = parse_sequence_short(p, n, seq_can, &chain->seq, verbose2);
	// read the canonical sequence 
	if ( esl_memstrcmp(p, n, scomma)  && readnow) parse_sequence_end (p, n, scomma, &readnow,    chain,     verbose2);  // sequence ends with ';'
	if (!esl_memstrpfx(p, n, seq_can) && readnow) parse_sequence_long(p, n,         &chain->seq, &chain->L, FALSE);
	if (esl_memstrpfx(p, n, seq_chainname)) {
	  // the sequence finished. Add as a chain to the pdbx structure.
	  // If different chains use the same structure, we will add more chains later.
	  rview_pdbx_AddChain(pdbx, chain, -1, verbose);
	  rview_chain_Destroy(1, chain); chain = NULL;
	}

      }
      else {  // if multiseq grab everything after 'targetid' and before the next '#' to get the chains
	if ( esl_memstrpfx(p, n, seq_targetid))          readnow = TRUE;
	if ( esl_memstrpfx(p, n, endloop)    && readnow) parse_add_allseqs(pdbx, &readnow, chunk, len, FALSE); // have to place between the other two options
	if (!esl_memstrpfx(p, n, newcommand) && readnow) parse_multiseq   (p, n, &chunk,         &len, FALSE);
      }

      // check if any of the sequences has heterogeneous residues.
      // That is positions with more than one residue residues.
      // This alternative residues have different atom coordinates
      // For heterogeneous residues, I consider all atoms
      if      ( esl_memstrpfx(p, n, seq_hetero))              checkhetero = TRUE;
      else if ( esl_memstrpfx(p, n, endloop)  && checkhetero) checkhetero = FALSE;
      else if (checkhetero) parse_sequence_hetero(pdbx, p, n, TRUE);
		      
      // Figure out how many chains per sequence and their names.
      //
      // The chain names _struct_asym.id appear in a loop.
      // field _struct_asym.entity_id identifies the sequence corresponding to that chain.
      //
      // Later on, in ATOMS, the field atom_site.label_asym_id is the name of the chain found in _struc_asym.id
      // _atom_site.label_asym_id DOESN'T correspond with _entity_poly.pdbx_strand_id
      //     (for _atom_site.label_entity_id == _entity_poly.id) 
      //
      //
      // #
      // loop_
      // _struct_asym.id 
      // _struct_asym.pdbx_blank_PDB_chainid_flag 
      // _struct_asym.pdbx_modified 
      // _struct_asym.entity_id 
      // _struct_asym.details 
      // A N N 1 ? 
      // B N N 2 ? 
      // C N N 2 ? 
      //
      // which tells us there is one chain (A) for sequence 1 and two chains (B,C) for sequence 2.
      //
      if ( esl_memstrpfx(p, n, sasym))              { parse_sasymfield(      p, n, nsl, sasymlabel, &nsaf, &usesasymf, FALSE); issasym = TRUE; }
      if ( esl_memstrpfx(p, n, endloop) && issasym) issasym = FALSE; // have to place in between the two other options
      if (!esl_memstrpfx(p, n, sasym)   && issasym) { parse_sasym     (pdbx, p, n,                  nsaf,  usesasymf,  verbose2);              }
      
      // THE ATOMS
      //
      // we are interested in 8 of the fields, which can appear in any ordering.
      // First, we need to figure out the ordering for this particular file.
      //
      // # 
      // loop_
      // _atom_site.group_PDB             (type)
      // _atom_site.id                    (idx)
      // _atom_site.type_symbol           (atomid)
      // _atom_site.label_atom_id         
      // _atom_site.label_alt_id 
      // _atom_site.label_comp_id         (reschr)
      // _atom_site.label_asym_id         (chain)
      // _atom_site.label_entity_id       (seqid)
      // _atom_site.label_seq_id          (resid)
      // _atom_site.pdbx_PDB_ins_code 
      // _atom_site.Cartn_x               (x)
      // _atom_site.Cartn_y               (y)
      // _atom_site.Cartn_z               (z)
      // _atom_site.occupancy 
      // _atom_site.B_iso_or_equiv 
      // _atom_site.pdbx_formal_charge 
      // _atom_site.auth_seq_id 
      // _atom_site.auth_comp_id 
      // _atom_site.auth_asym_id 
      // _atom_site.auth_atom_id 
      // _atom_site.pdbx_PDB_model_num
      //
      // examples,
      //
      // HETATM 1    P  PG    . GTP A 1 1  ? 12.672  45.003 115.432 0.50 18.30  ? 6   GTP A PG    1 
      // ATOM   33   P  P     . G   A 1 2  ? 7.252   42.695 107.844 1.00 43.04  ? 7   G   A P     1 
      //
      if (esl_memstrpfx(p, n, atomsite)) parse_atomfield(      p, n, nal, atmlabel, &natf, &useatmf, FALSE);
      if (esl_memstrpfx(p, n, atomid))   parse_atom     (pdbx, p, n, atomid,        natf,  useatmf,  verbose);
      
      // HETATM are defined as not beeing attached to anything, so I would not include them
      //if (esl_memstrpfx(p, n, hatomid))  parse_atom     (pdbx, p, n, hatomid, natf, useatmf, TRUE);
    }
 
  if (status != eslEOF) esl_fatal("file %s: expected EOF, got code %d", bf->filename, status);
  esl_buffer_Close(bf);

  // we got it. Print if asked
  if (1||verbose) rview_pdbx_Write(stdout, pdbx, FALSE);

  // Check if the res numbering agrees with the location in the sequence
  //
  // in .pbd format that is often no the case, and a pain in the neck.
  // in .cif format it appears to be better.
  //
  // If the pdb file does not checksum, we will make a smith-waterman alignment
  // a create a map of the sequence to the resseq.
  if (rview_pdbx_Checksum(pdbx, verbose) != eslOK) esl_fail("pdbx did not checksum");
  
  // create the map between chain->seq and the sequence described by the atoms
  //
  
  *ret_pdbx = pdbx;

  // clean up
  free(chunk);
  free(usesasymf);
  free(useatmf);
  for (s = 0; s < nsl; s ++) free(sasymlabel[s]);
  free(sasymlabel);
  for (l = 0; l < nal; l ++) free(atmlabel[l]);
  free(atmlabel);
  if (chain) rview_chain_Destroy(1, chain);
  return eslOK;
  
 ERROR:
  if (chunk)     free(chunk);
  if (usesasymf) free(usesasymf);
  if (useatmf)   free(useatmf);
  for (s = 0; s < nsl; s ++) if (sasymlabel[s]) free(sasymlabel[s]);
  if (sasymlabel) free(sasymlabel);
  for (l = 0; l < nal; l ++) if (atmlabel[l]) free(atmlabel[l]);
  if (atmlabel) free(atmlabel);
  if (chain) rview_chain_Destroy(1, chain);
  if (pdbx)  rview_pdbx_Destroy(pdbx);
  return status;
}

void
rview_atom_Init(ATOM *atom)
{
  atom->type     = ATYPE_NONE;
  atom->idx      = -1;
  atom->atomid   = NULL;
  atom->atomidx  = NULL;
  atom->chain    = NULL;
  atom->reschr   = NULL;
  atom->seqid    = -1;
  atom->resid    = -1;
  atom->residx   = -1;
}

void
rview_atom_Destroy(esl_pos_t na, ATOM *atom)
{
  ATOM      *each;
  esl_pos_t  a;
  
  if (atom) {
    for (a = 0; a < na; a ++) {
      each = &atom[a];
      if (each->atomid)  free(each->atomid);
      if (each->atomidx) free(each->atomidx);
      if (each->chain)   free(each->chain);
      if (each->reschr)  free(each->reschr);
    }
    free(atom);
  }
}

void
rview_atom_Write(FILE *fp, ATOM *atom)
{
  if (atom == NULL) return;
  fprintf(fp, "%lld %s %lld %s %s %lld %s %s %f %f %f\n",
	  atom->seqid, atom->chain, atom->resid, atom->reschr, (atom->type == ATYPE_ATOM)? "ATOM":"HETATM", atom->idx,
	  atom->atomid, atom->atomidx, atom->x, atom->y, atom->z);
}

struct chain_s *
rview_chain_Create()
{
  struct chain_s *chain = NULL;
  int             status;

  ESL_ALLOC(chain, sizeof(struct chain_s));
  chain->name      = NULL;
  chain->seqid     = -1;
  chain->L         = 0;
  chain->hetero    = FALSE;
  chain->seq       = NULL;
  chain->seqtype   = NULL;
  chain->nr        = 0;
  chain->res       = NULL;
  chain->resseq    = NULL;
  chain->atom_map  = NULL;
  chain->res       = NULL;
  
  return chain;

 ERROR:
  return NULL;
}

int
rview_chain_AddRes(struct chain_s *chain)
{
  int status;
  if (chain->nr == 0) ESL_ALLOC  (chain->res, sizeof(RES));
  else                ESL_REALLOC(chain->res, sizeof(RES) * (chain->nr+1));
  rview_res_Init(&chain->res[chain->nr]);

  chain->nr ++;
  return eslOK;

 ERROR:
  return eslFAIL;
}

void
rview_chain_Destroy(esl_pos_t nc, struct chain_s *chain)
{
  struct chain_s *each;
  esl_pos_t       c;
  
  if (chain) {
    for (c = 0; c < nc; c ++) {
      if (each == NULL) continue;
      each = &chain[c];
      if (each->name)    free(each->name);
      if (each->seqtype) free(each->seqtype);
      
      if (each->seq)    free(each->seq);
      if (each->resseq) free(each->resseq);
      
      rview_res_Destroy(each->nr, each->res);
      if (each->atom_map) free(each->atom_map);
    }
    free(chain);
  }
}

void
rview_chain_Write(FILE *fp, struct chain_s *chain, int verbose)
{
  int r;

  if (chain == NULL) return;
  
  if (chain->name)    fprintf(fp, "\nchain: %s\n", chain->name); 
  if (chain->seqtype) fprintf(fp, "type: %s\n",    chain->seqtype); 
  if (chain->hetero)  fprintf(fp, "seqid: %lld (hetero)\n", chain->seqid);
  else                fprintf(fp, "seqid: %lld\n", chain->seqid); 
  if (chain->seq)     fprintf(fp, "sequence: (%lld)\n%s\n", chain->L, chain->seq);

  if (chain->resseq)
    fprintf(fp, "res_sequence: (%lld)\n%s\n", chain->nr, chain->resseq);
  else
    fprintf(fp, "nres: (%lld)\n", chain->nr);
  
  for (r = 0; r < chain->nr; r ++)
    rview_res_Write(fp, &chain->res[r], verbose);
}

int
rview_pdbx_AddAtom(PDBX *pdbx, ATOM *atom, int verbose)
{
  struct chain_s *chain;
  RES            *res;
  int             c;
  int             r = 0;
  int             status;
  
  if (atom == NULL) return eslOK;

  // find the chain
  // chains have already been allocated and filled with the sequences
  for (c = 0; c < pdbx->nch; c++) {
    chain = &(pdbx->chain[c]);
    if (esl_strcmp(atom->chain, chain->name) == eslOK) break;
  }
  if (c == pdbx->nch) esl_fatal("could not find chain %s in pdbx %s\n", atom->chain, pdbx->pdbname);

  // find the residue or create a new one
  for (r = 0; r < chain->nr; r ++)  
    if (atom->resid == chain->res[r].resnum) break;

  if (chain->nr == 0 || r == chain->nr) {
    rview_chain_AddRes(chain);
    res = &(chain->res[r]);
    res->na           = 0;
    res->h            = 1;
    res->seqid        = atom->seqid;
    res->resnum       = atom->resid;
    res->from_atomidx = atom->idx;
    ESL_ALLOC(res->reschr, sizeof(char *));
    esl_sprintf(&res->reschr[0], atom->reschr);
    esl_strcat(&chain->resseq, -1, atom->reschr, -1); // add the residue to chain->resseq
    esl_sprintf(&res->chain,  atom->chain);
  }

  // now add the atom to this residue
  rview_res_AddAtom(&(chain->res[r]), atom);
  if (chain->res[r].h > 1) chain->hetero = TRUE;

  return eslOK;

 ERROR:
  return status;
}

int
rview_pdbx_AddChain(PDBX *pdbx, struct chain_s *chain, int idx, int verbose)
{
  struct chain_s  *new;
  int              c;
  int              status;

  if (chain == NULL) return eslOK;
  if (idx > pdbx->nch) idx = pdbx->nch; // attach at the end
  if (idx < 0)         idx = pdbx->nch; // attach at the end

  if (pdbx->chain == NULL) ESL_ALLOC  (pdbx->chain, sizeof(struct chain_s) * (pdbx->nch+1));
  else                     ESL_REALLOC(pdbx->chain, sizeof(struct chain_s) * (pdbx->nch+1));

  for (c = pdbx->nch; c > idx; c --) pdbx->chain[c] = pdbx->chain[c-1];
  
  new           = &(pdbx->chain[idx]);
  new->L        = strlen(chain->seq);
  new->hetero   = chain->hetero;
  new->seqid    = chain->seqid;
  new->nr       = 0;
  new->res      = NULL;
  new->resseq   = NULL;
  new->atom_map = NULL;
  
  esl_sprintf(&new->name,    chain->name);
  esl_sprintf(&new->seqtype, chain->seqtype);
  esl_sprintf(&new->seq,     chain->seq);
  
  pdbx->nch ++;
  if (verbose) rview_chain_Write(stdout, new, verbose);
  
   return eslOK;
  
 ERROR:
  return status;
}

// This function checks for all chains whether
// the resnum assigned in the atoms, corresponds to
// the residue's position in the sequence.
//
// on success, it returns eslOK, otherwise eslFAIL;
// 
int
rview_pdbx_Checksum(PDBX *pdbx, int verbose)
{
  struct chain_s *chain;
  RES            *res;
  char           *seq;
  char            seqchr[2];
  int             c;
  int             r;
  int             h;
  int             checksum = eslFAIL;
  
  for (c = 0; c < pdbx->nch; c ++) {
    chain = &(pdbx->chain[c]);
    seq   = chain->seq;
    
    for (r = 0; r < chain->nr; r ++) {
      res = &(chain->res[r]);
      
      seqchr[0] = seq[res->resnum-1];
      seqchr[1] = '\0';
      for (h = 0; h < res->h; h ++) 
	checksum = esl_strcmp(res->reschr[h], seqchr);
      if (checksum != eslOK) {
	if (verbose)
	  printf("chain %s (seqid %lld) does not checksum at residue %lld (seqres %s)\n",
		 chain->name, chain->seqid, res->resnum, seqchr); 
      }
    }
    
  }
  
  return checksum;
}

PDBX *
rview_pdbx_Create()
{
  PDBX *pdbx = NULL;
  int   status;

  ESL_ALLOC(pdbx, sizeof(PDBX));
  
  pdbx->pdbname  = NULL;
  pdbx->nsq      = 0;
  pdbx->nch      = 0;
  pdbx->chain    = NULL;
  pdbx->res_low  = 0.;
  pdbx->res_high = 0.;
  
  return pdbx;
  
 ERROR:
  return NULL;
}

void
rview_pdbx_Destroy(PDBX *pdbx)
{
  if (pdbx) {
    if (pdbx->pdbname) free(pdbx->pdbname);
    if (pdbx->chain)   rview_chain_Destroy(pdbx->nch, pdbx->chain);
    free(pdbx);
  }
}

void
rview_pdbx_Write(FILE *fp, PDBX *pdbx, int verbose)
{
  int c;
  fprintf(fp, "\nPDB:      %s\n", pdbx->pdbname);
  fprintf(fp, "resolution: %f %f\n", pdbx->res_low, pdbx->res_high); 
  fprintf(fp, "nsequences: %d\n", pdbx->nsq);
  fprintf(fp, "nchains:    %d\n", pdbx->nch);
  for (c = 0; c < pdbx->nch; c ++)
    rview_chain_Write(fp, &pdbx->chain[c], verbose);
}

int
rview_res_AddAtom(RES *res, ATOM *atom)
{
  ATOM *new;
  int   h;
  int   found = FALSE;
  int   status;

  if (res->na == 0) ESL_ALLOC  (res->atom, sizeof(ATOM));
  else              ESL_REALLOC(res->atom, sizeof(ATOM) * (res->na+1));

  if (res->na > 0) { //check the atom belongs to this residue
    if (atom->seqid != res->seqid)                      esl_fatal("atom does not belong to this residue. bad seqid");
    if (atom->resid != res->resnum)                     esl_fatal("atom does not belong to this residue. bad resid");
    if (esl_strcmp(atom->chain,  res->chain)  != eslOK) esl_fatal("atom does not belong to this residue. bad chain");
    
    // the residue agrees? (don't forget heterogeneity)
    for (h = 0; h < res->h; h++) 
      if (esl_strcmp(atom->reschr, res->reschr[h]) == eslOK) { found = TRUE; break; }
    
    if (found == FALSE) { // heterogeneous residue
      ESL_REALLOC(res->reschr, sizeof(char *) * (res->h+1));
      esl_sprintf(&res->reschr[res->h], atom->reschr);
      res->h ++;
      printf("hetero residue %lld in chain %s seq %lld\n", res->resnum, atom->chain, atom->seqid);
    }

  }
  
  new = &(res->atom[res->na]);
  new->type   = atom->type;
  new->idx    = atom->idx;
  new->seqid  = atom->seqid;
  new->resid  = atom->resid;
  esl_sprintf(&new->atomid,  atom->atomid);
  esl_sprintf(&new->atomidx, atom->atomidx);
  esl_sprintf(&new->chain,   atom->chain);
  esl_sprintf(&new->reschr,  atom->reschr);
  
  res->na ++;
  return eslOK;
  
 ERROR:
  return eslFAIL;
}

void
rview_res_Init(RES *res)
{
  res->h      = 1;
  res->reschr = NULL;
  res->chain  = NULL;

  res->from_atomidx = -1;
  res->resnum       = -1;
  
  res->na   = 0;
  res->atom = NULL;
}

void
rview_res_Destroy(esl_pos_t nr, RES *res)
{
  RES       *each;
  esl_pos_t  r;
  esl_pos_t  h;

  if (res) {
    for (r = 0; r < nr; r ++) {
      each = &res[r];
      for (h = 0; h < each->h; h ++) if (each->reschr[h]) free(each->reschr[h]);
      if (each->reschr) free(each->reschr);
      if (each->chain)  free(each->chain);
      if (each->atom)   rview_atom_Destroy(each->na, each->atom);
    }
    free(res);
  }
}

void
rview_res_Write(FILE *fp, RES *res, int verbose)
{
  int a;

  if (res == NULL) return;
  
  fprintf(fp, "%s %lld %s natoms %lld idx %lld\n", res->reschr[0], res->resnum, res->chain, res->na, res->from_atomidx);
  if (verbose) {
    for (a = 0; a < res->na; a ++)
      rview_atom_Write(fp, &res->atom[a]);
  }
}


/*------------------------------------- static routines ---------------------------------------*/


// The struct_asym
//
// gives use the assotiation between the chains () and the sequences ()

// we are interested in 2 'labels' of all possible fields 
//
// There is not a unique order of all the fields,
// we need to figure out in which order they appeared. 

// # 
// loop_

// #
//
//
//

static int
sasymlabel_create(int *ret_nsl, char ***ret_sasymlabel, int verbose)
{
  char **sasymlabel = NULL;
  char   sasym_id[]       = "_struct_asym.id";
  char   sasym_entityid[] = "_struct_asym.entity_id";
  int    nsl = 2;
  int    s = 0;
  int    status;

  // We want to keep 2 of all possible _struct_asym labels
  //
  ESL_ALLOC(sasymlabel, sizeof(char *) * nsl);
  esl_sprintf(&sasymlabel[s++], sasym_id);
  esl_sprintf(&sasymlabel[s++], sasym_entityid);

  if (verbose) 
    for (s = 0; s < nsl; s ++) printf("struct_asym_label %s\n", sasymlabel[s]);
  
  *ret_nsl        = nsl;
  *ret_sasymlabel = sasymlabel;
  return eslOK;
  
 ERROR:
  for (s = 0; s < nsl; s ++)
    if (sasymlabel[s]) free(sasymlabel[s]);
  if (sasymlabel) free(sasymlabel);
  return status;
}

// THE ATOMS
//
// we are interested in 8 'labels' of all possible fields 
//
// There is not a unique order of all the fields,
// we need to figure out in which order they appeared. 

// # 
// loop_
// _atom_site.group_PDB             (type,   0)
// _atom_site.id                    (idx,    1)
// _atom_site.type_symbol           (atmid,  2)
// _atom_site.label_atom_id         (atmidx  3)
// _atom_site.label_alt_id 
// _atom_site.label_comp_id         (reschr, 5)
// _atom_site.label_asym_id         (chain,  6)
// _atom_site.label_entity_id       (seqid,  7)
// _atom_site.label_seq_id          (resid,  8      
// _atom_site.pdbx_PDB_ins_code 
// _atom_site.Cartn_x               (x,      10)
// _atom_site.Cartn_y               (y,      11)
// _atom_site.Cartn_z               (z,      12)
// _atom_site.occupancy 
// _atom_site.B_iso_or_equiv 
// _atom_site.pdbx_formal_charge 
// _atom_site.auth_seq_id 
// _atom_site.auth_comp_id 
// _atom_site.auth_asym_id 
// _atom_site.auth_atom_id 
// _atom_site.pdbx_PDB_model_num
// #
//
// examples,
//
// HETATM 1    P  PG    . GTP A 1 1  ? 12.672  45.003 115.432 0.50 18.30  ? 6   GTP A PG    1 
// ATOM   33   P  P     . G   A 1 2  ? 7.252   42.695 107.844 1.00 43.04  ? 7   G   A P     1 
//
// I create a map
//
//
//
static int
atmlabel_create(int *ret_nal, char ***ret_atmlabel, int verbose)
{
  char **atmlabel = NULL;
  char   atom_type[]   = "_atom_site.group_PDB";
  char   atom_idx[]    = "_atom_site.id";
  char   atom_chr[]    = "_atom_site.type_symbol";
  char   atom_chrx[]   = "_atom_site.label_atom_id";
  char   atom_reschr[] = "_atom_site.label_comp_id";
  char   atom_chain[]  = "_atom_site.label_asym_id";
  char   atom_seqid[]  = "_atom_site.label_entity_id";
  char   atom_resid[]  = "_atom_site.label_seq_id";
  char   atom_xcoord[] = "_atom_site.Cartn_x";
  char   atom_ycoord[] = "_atom_site.Cartn_y";
  char   atom_zcoord[] = "_atom_site.Cartn_z";
  int    nal = 11;
  int    l = 0;
  int    status;

  // this are the 9 fields we want to keep
  ESL_ALLOC(atmlabel, sizeof(char *) * nal);
  esl_sprintf(&atmlabel[l++], atom_type);
  esl_sprintf(&atmlabel[l++], atom_idx);
  esl_sprintf(&atmlabel[l++], atom_chr);
  esl_sprintf(&atmlabel[l++], atom_chrx);
  esl_sprintf(&atmlabel[l++], atom_reschr);
  esl_sprintf(&atmlabel[l++], atom_chain);
  esl_sprintf(&atmlabel[l++], atom_seqid);
  esl_sprintf(&atmlabel[l++], atom_resid);
  esl_sprintf(&atmlabel[l++], atom_xcoord);
  esl_sprintf(&atmlabel[l++], atom_ycoord);
  esl_sprintf(&atmlabel[l++], atom_zcoord);

  if (verbose) 
    for (l = 0; l < nal; l ++) printf("atmlabel %s\n", atmlabel[l]);
  
  *ret_nal       = nal;
  *ret_atmlabel = atmlabel;
  return eslOK;
  
 ERROR:
  for (l = 0; l < nal; l ++)
    if (atmlabel[l]) free(atmlabel[l]);
  if (atmlabel) free(atmlabel);
  return status;
}

static int
parse_atomfield(char *p, esl_pos_t n, int nal, char **atmlabel, int *ret_natf, int **ret_useatmf, int verbose)
{
  int  natf    = *ret_natf;
  int *useatmf = *ret_useatmf;
  int  l;
  int  status;

  if (natf == 0) ESL_ALLOC  (useatmf, sizeof(int));
  else            ESL_REALLOC(useatmf, sizeof(int) * (natf+1));
  
  useatmf[natf] = -1; // atom field not used by default
  for (l = 0; l < nal; l ++) {
    if  (esl_memstrpfx(p, n, atmlabel[l])) {
      useatmf[natf] = l;
      if (verbose) printf ("%s is atmfield %d\n", atmlabel[l], natf);
    }
  }
  
  *ret_natf    = ++natf;
  *ret_useatmf = useatmf;
  return eslOK;

 ERROR:
  if (useatmf) free(useatmf);
  return status;
}

// ATOM   33   P  P     . G   A 1 2  ? 7.252   42.695 107.844 1.00 43.04  ? 7   G   A P     1 
// HETATM 1    P  PG    . GTP A 1 1  ? 12.672  45.003 115.432 0.50 18.30  ? 6   GTP A PG    1 
static int
parse_atom(PDBX *pdbx, char *p, esl_pos_t n, char *atomid, int naf, int *useatmfield, int verbose)
{
  ATOM        *atom = NULL;
  char      **field = NULL;
  esl_pos_t   i;
  esl_pos_t   len = 0;
  esl_pos_t   salloc = n + 1;
  int         idx = 0;
  int         isblank = FALSE;
  int         f;
  int         status;

  if (! esl_memstrpfx(p, n, atomid)) esl_fail("not an atom %s\n", atomid);
  
  ESL_ALLOC(atom,     sizeof(ATOM));
  ESL_ALLOC(field,    sizeof(char *) * naf);
  ESL_ALLOC(field[0], sizeof(char  ) * naf * salloc);
  for (f = 1; f < naf; f ++) field[f] = field[f-1] + salloc;
  
  for (i = 0; i < n; i++) {
    if (! isspace(p[i])) {
      isblank = FALSE;
      field[idx][len++] = p[i];
    }
    else {
      if (isblank == FALSE) { field[idx][len] = '\0';  idx ++; }
      isblank = TRUE;
      len     = 0;
    }
  }
  if (idx != naf) esl_fatal("wrong number of atom fields %d should be %d\n", idx, naf);  

  // asign the fields to the ATOM structure (hard coded)
  for (f = 0; f < naf; f++) {
    if (useatmfield[f] == 0) {
      if      (esl_strcmp(field[f],"ATOM")   == eslOK) atom->type = ATYPE_ATOM;
      else if (esl_strcmp(field[f],"HETATM") == eslOK) atom->type = ATYPE_HETATM;
      else esl_fatal("wrong type of atom %s\n", field[f]);
    }
    
    if (useatmfield[f] ==  1) atom->idx    = atoi(field[f]);
    if (useatmfield[f] ==  2) esl_sprintf  (&atom->atomid,  field[f]);
    if (useatmfield[f] ==  3) esl_sprintf  (&atom->atomidx, field[f]);
    if (useatmfield[f] ==  4) res_oneletter(&atom->reschr,  field[f]);
    if (useatmfield[f] ==  5) esl_sprintf  (&atom->chain,   field[f]);
    
    if (useatmfield[f] ==  6) atom->seqid = atoi(field[f]);
    if (useatmfield[f] ==  7) atom->resid = atoi(field[f]);
    if (useatmfield[f] ==  8) atom->x     = atof(field[f]);
    if (useatmfield[f] ==  9) atom->y     = atof(field[f]);
    if (useatmfield[f] == 10) atom->z     = atof(field[f]);
  }
  
  if (verbose) rview_atom_Write(stdout, atom);
  
  rview_pdbx_AddAtom(pdbx, atom, verbose);

  rview_atom_Destroy(1, atom);
  free(field[0]);
  free(field);
  return eslOK;

 ERROR:
  if(field[0]) free(field[0]);
  if (field)   free(field);
  if (atom)    rview_atom_Destroy(1, atom);
  return status;
}

static int
parse_pdbname(char *p, esl_pos_t n, char *entryid, char **ret_pdbname, int verbose)
{
  char       *pdbname = NULL;
  esl_pos_t   slen;
  esl_pos_t   len = 0;
  esl_pos_t   i;
  int         status;
  
  if (! esl_memstrpfx(p, n, entryid)) esl_fail("not an entryid %s\n", entryid);
  
  slen = (esl_pos_t)strlen(entryid);
  ESL_ALLOC(pdbname, sizeof(char)*(n-slen+1));
  
  for (i = slen; i < n; i++) 
    if (! isspace(p[i])) pdbname[len++] = p[i]; 
  pdbname[len] = '\0';
  
  if (verbose) printf("pdbname %s\n", pdbname);
  
  *ret_pdbname = pdbname;
  
  return eslOK;

  ERROR:
  if (pdbname) free(pdbname);
  return status;
}

static
int parse_reslow(char *p, esl_pos_t n, char *rlow, double *ret_res_low, int verbose)
{
  char       *resolution = NULL;
  esl_pos_t   slen;
  esl_pos_t   salloc;
  esl_pos_t   len = 0;
  esl_pos_t   i;
  int         status;
  
 if (! esl_memstrpfx(p, n, rlow)) esl_fail("not a resolution_low entry %s\n", rlow);
     
  slen = (esl_pos_t)strlen(rlow);
  salloc = n - slen + 1;
  
  if (salloc > 0) {
    ESL_ALLOC(resolution, sizeof(char) * salloc);
    for (i = slen; i < n; i++) 
      if (! isspace(p[i])) resolution[len++] = p[i]; 
    resolution[len] = '\0';

    if (ret_res_low) {
      *ret_res_low = atof(resolution);
      if (verbose) printf("res_low %f\n", *ret_res_low);
    }
    
    free(resolution);
  }
  return eslOK;

 ERROR:
  if (resolution) free(resolution);
  return status;
}

static
int parse_reshigh(char *p, esl_pos_t n, char *rhigh, double *ret_res_high, int verbose)
{
  char       *resolution = NULL;
  esl_pos_t   slen;
  esl_pos_t   salloc;
  esl_pos_t   len = 0;
  esl_pos_t   i;
  int         status;
  
  if (! esl_memstrpfx(p, n, rhigh)) esl_fail("not a resolution_high entry %s\n", rhigh);
  
  slen   = (esl_pos_t)strlen(rhigh);
  salloc = n - slen + 1;

  if (salloc > 0) {
    ESL_ALLOC(resolution, sizeof(char) * salloc);
    for (i = slen; i < n; i++) 
      if (! isspace(p[i])) resolution[len++] = p[i]; 
    resolution[len] = '\0';
    
    *ret_res_high = atof(resolution);
    if (verbose) printf("res_high %f\n", *ret_res_high);
    
    free(resolution);
  }
  return eslOK;

  ERROR:
  if (resolution) free(resolution);
  return status;

}

static
int parse_seqtype(char *p, esl_pos_t n, char *type, char **ret_seqtype, int verbose)
{
  char       *seqtype = NULL;
  esl_pos_t   slen;
  esl_pos_t   len = 0;
  esl_pos_t   i;
  int         status;
  
  if (! esl_memstrpfx(p, n, type)) esl_fail("not an seqtype %s\n", type);
  
  slen = (esl_pos_t)strlen(type);
  ESL_ALLOC(seqtype, sizeof(char)*(n-slen+1));
  
  for (i = slen; i < n; i++) 
    if (! isspace(p[i])) seqtype[len++] = p[i]; 
  seqtype[len] = '\0';
  
  if (verbose) printf("seqtype %s\n", seqtype);
  
  *ret_seqtype = seqtype;
  
  return eslOK;

  ERROR:
  if (seqtype) free(seqtype);
  return status;
}

static
int parse_seqid(PDBX *pdbx, char *p, esl_pos_t n, char *entity, int64_t *ret_seqid, int *ret_multiseq, int verbose) {
  char       *seqid = NULL;
  int         multiseq = TRUE;
  esl_pos_t   slen;
  esl_pos_t   len = 0;
  esl_pos_t   salloc = 0;
  esl_pos_t   i;
  int         status;
  
  if (! esl_memstrpfx(p, n, entity)) esl_fail("not an entity %s\n", entity);
  
  slen = (esl_pos_t)strlen(entity);
  for (i = slen; i < n; i++) 
    if (! isspace(p[i])) salloc ++;

  // single sequence
  if (salloc > 0) {
    multiseq = FALSE;
    
    ESL_ALLOC(seqid, sizeof(char) * (salloc+1));
    for (i = slen; i < n; i++) 
      if (! isspace(p[i])) {
	seqid[len++] = p[i];
      }
    seqid[len] = '\0';
    pdbx->nsq ++;
  }
  if (verbose) printf("multiseq? %d\n", multiseq);

  *ret_seqid      = (seqid)? (int64_t)atoi(seqid) : -1;
  *ret_multiseq = multiseq;

  free(seqid);
  return eslOK;

 ERROR:
  if (seqid) free(seqid);
  return status;
}

static int
parse_sequence_short(char *p, esl_pos_t n, char *seq_can, char **ret_seq, int verbose)
{
  int         readnow = FALSE;
  char       *seq     = NULL;
  esl_pos_t   slen;
  esl_pos_t   len = 0;
  esl_pos_t   salloc;
  esl_pos_t   i;
  int         status;
  
  if (! esl_memstrpfx(p, n, seq_can)) esl_fail("not a seq_can %s\n", seq_can);

  slen = (esl_pos_t)strlen(seq_can);
  salloc = 0;
  for (i = slen; i < n; i++) 
    if (! isspace(p[i])) salloc ++;

  if (salloc == 0) readnow = TRUE;
  else {
    // short sequence fits in this line
    ESL_ALLOC  (seq, sizeof(char) * salloc);
    for (i = slen; i < n; i++) 
      if (! isspace(p[i])) seq[len++] = p[i];
    seq[len] = '\0';
  }

  *ret_seq = seq;
  return readnow;

 ERROR:
  *ret_seq = NULL;
  return readnow;
}
static
int parse_sequence_long(char *p, esl_pos_t n, char **ret_seq, esl_pos_t *ret_seqlen, int verbose)
{
  char       *seq    = *ret_seq;
  esl_pos_t   seqlen = *ret_seqlen;
  esl_pos_t   salloc;
  esl_pos_t   i;
  int         status;

  salloc = seqlen + n + 1;  
  if (seq == NULL) ESL_ALLOC  (seq, sizeof(char) * salloc);
  else             ESL_REALLOC(seq, sizeof(char) * salloc);  
  
  for (i = (p[0]==';')? 1:0; i < n; i++) {
    if (! isspace(p[i])) { seq[seqlen++] = p[i]; }
  }

  if (verbose) printf("canonical sequence\n%s\n", seq);

  *ret_seqlen = seqlen;
  *ret_seq    = seq;
  return eslOK;

 ERROR:
  if (seq) free(seq);
  return status;
}

static int
parse_sequence_hetero(PDBX *pdbx, char *p, esl_pos_t n, int verbose)
{
  char      **field = NULL;
  char       *reschr;
  esl_pos_t   seqid;
  esl_pos_t   resnum;
  esl_pos_t   i;
  esl_pos_t   len = 0;
  esl_pos_t   salloc = n + 1;
  int         idx = 0;
  int         isblank = FALSE;
  esl_pos_t   nf;
  esl_pos_t   f;
  esl_pos_t   c, r;
  int         status;

  ESL_ALLOC(field,      sizeof(char *));
  ESL_ALLOC(field[idx], sizeof(char) * salloc);
  for (i = 0; i < n; i++) {
    if (! isspace(p[i])) {
      isblank = FALSE;
      field[idx][len++] = p[i];
    }
    else {
      if (isblank == FALSE) {
	field[idx][len] = '\0';
	idx ++;
	ESL_REALLOC(field,      sizeof(char *) * (idx+1));
	ESL_ALLOC  (field[idx], sizeof(char)   * salloc);
      }
      isblank = TRUE;
      len     = 0;
    }
  }
  nf = idx;
  if (nf != 4) esl_fatal("sequence_hetero should have 4 fields not %d\n", nf);

  seqid  = atoi(field[0]);
  resnum = atoi(field[1]);
  reschr = field[2];
  if (esl_strcmp(field[3], "y") == eslOK) {
    if (verbose) printf("hetero sequence %lld residue %lld chr %s\n", seqid, resnum, reschr);    
    for (c = 0; c < pdbx->nch; c ++) 
      if (pdbx->chain[c].seqid == seqid)
	pdbx->chain[c].hetero = TRUE;     
  }
  
  for (f = 0; f <= nf; f ++) free(field[f]);
  free(field);
  return eslOK;

 ERROR:
  for (f = 0; f <= nf; f ++) if (field[f]) free(field[f]);
  if (field) free(field);
  return status;
}

static
int parse_sequence_end(char *p, esl_pos_t n, char *scomma, int *ret_readnow, struct chain_s *chain, int verbose)
{
  int  readnow = FALSE;

  if (! esl_memstrpfx(p, n, scomma)) esl_fail("not a sequence end %s\n", scomma);
  
  chain->seq[chain->L] = '\0';
  if (verbose) printf("canonical sequence\n%s\n", chain->seq);

  *ret_readnow = readnow;
  
  return eslOK;
}


     
static int
parse_multiseq(char *p, esl_pos_t n, char **ret_chunk, esl_pos_t *ret_len, int verbose)
{
  char       *chunk = *ret_chunk;
  esl_pos_t   len    = *ret_len;
  esl_pos_t   salloc = len+n+2;
  esl_pos_t   i;
  int         status;

  if (chunk == NULL) ESL_ALLOC  (chunk, sizeof(char) * salloc);
  else               ESL_REALLOC(chunk, sizeof(char) * salloc);  

  for (i = 0; i < n; i++) chunk[len++] = p[i];
  chunk[len++] = '\n';
  
  if (verbose) printf("%s\n", chunk);

  *ret_len   = len;
  *ret_chunk = chunk;
  return eslOK;

 ERROR:
  if (chunk) free(chunk);
  return status;
}

/* pdbx provides sequences. Each sequence corresponds to one or more chains
 * 
 */
static int
parse_add_allseqs(PDBX *pdbx, int *ret_readnow, char *chunk, esl_pos_t len, int verbose)
{
  struct chain_s *chain = NULL;
  char           *tok;
  SEQSTATUS       seqstatus = SEQ_NONE;
  esl_pos_t       i;
  int             isshort;
  int             readnow  = FALSE;
  int             verbose2 = TRUE;
  int             status;

  chunk[len] = '\0'; // termine chunk as a proper string

  chain = rview_chain_Create();
  
  while ( (status = esl_strtok(&chunk, "\n", &tok)) == eslOK)
    {
      // a new sequence starts with an integer which is the seq_id
      // but watchout a chain can also be a number! and be at the begining of a line
      if (chunk_is_newsequence(pdbx, tok, chain, &isshort, verbose2)) { // start a new chain
	seqstatus = SEQ_NONE;
	
	// sometimes if the sequence is very short,  all 8 fields fit in the first line
	// we  have stored the name, to know that the sequence has finished, but we don;t want to store the chain name
	//
	// 1  'polypeptide(L)'   no no MKTDFNQKIEQLKEFI MKTDFNQKIEQLKEFI A7    ? 
	//
	if (chain->name) {
	  free(chain->name); chain->name = NULL;
	  rview_pdbx_AddChain(pdbx, chain, -1, verbose);
	}
	continue;
      }
      
      // the short consensus sequence
      // format:
      //
      // sometimes is two lines breaking after seq:
      //
      //        chainnum seqtype no no seq
      //        seqcan chain
      //
      // example:
      //
      // 1  'polypeptide(L)'   no no MKTDFNQKIEQLKEFIEECRRVWLVLKKPTKDEYLAVAKVTALGISLLGIIGYIIHVPATYIKGILK 
      // MKTDFNQKIEQLKEFIEECRRVWLVLKKPTKDEYLAVAKVTALGISLLGIIGYIIHVPATYIKGILK A7    ? 
      //
      //
      // sometimes is two lines breaking after seqcan
      //
      //        chainnum seqtype no no seq seqcan
      //        chain
      //
      // example:
      //
      // 37 'polypeptide(L)'   no no MGQKWKLYEIKDGKVIRKNKFCPRCGPGVFMADHGDRWACGKCGYTEWKK MGQKWKLYEIKDGKVIRKNKFCPRCGPGVFMADHGDRWACGKCGYTEWKK 
      // BY    ? 
      //
      // Sometimes the sequence does not fit im one line , but the canonical sequence does
      // How is that possible, you may wonder? here is an example
      //
      // 85 'polypeptide(L)'   no no 
      // ;(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)
      // (UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)
      // (UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)(UNK)
      // ;
      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Ds    ? 
      //
      if (isshort) {
	chunk_is_sequence_short(tok, &chain->seq, NULL, verbose2);
	rview_pdbx_AddChain(pdbx, chain, -1, verbose);
      }
      else { 
	// the long consensus sequence
	// format:
	//
	// ;seq
	// seq_cont
	// ...
	// seq_cont
	// ;
	// ;seqcan
	// seqcan_cont
	// ...
	// seqcan_cont
	// ;
	//
	// or
	//
	// ;seq
	// seq_cont
	// ...
	// seq_cont
	// ;
	// seqcan AA ?
	// 
	// We want seqcan which has replaced all residues into canonical residues
	//
	if (chunk_is_sequence_start(tok, &chain->seq, &seqstatus)) {
	  if (verbose && seqstatus == SEQCAN_START) printf("seq_start %s\n", chain->seq);
	}
	if (chunk_is_sequence_cont (tok, &chain->seq, &seqstatus)) {
	  if (verbose && seqstatus == SEQCAN_CONT) printf("seq_cont %s\n", chain->seq);
	}
	if (chunk_is_sequence_end  (tok, &seqstatus)) {
	  if ((verbose2) && seqstatus == SEQCAN_END) printf("%s\n", chain->seq);
	}

	// After the sequence ends, we can add a chain.
	// AQ    ?
	//
	if (seqstatus == SEQCAN_END && chunk_is_chainname(tok, NULL)) {
	  // now we are ready to add the chain to pdbx.
	  // If several chains use the same sequence, we will add more chains with the same sequence later
	  isshort = FALSE;
	  rview_pdbx_AddChain(pdbx, chain, -1, verbose);
	}	
      }
    }
  if (status != eslEOL) esl_fatal("strtok failure");
  
  if (verbose2) printf("nseqs %d nch %d\n", pdbx->nsq, pdbx->nch);
  
  *ret_readnow = readnow;
  if (chain) rview_chain_Destroy(1, chain); 
  return eslOK; 
}

static int
chunk_is_newsequence(PDBX *pdbx, char *line, struct chain_s *chain, int *ret_isshort, int verbose)
{
  char      *tmp     = NULL;
  char      *seqtype = NULL;
  char      *tok;
  int        seqid;
  int        newseq = FALSE;
  int        isshort  = FALSE;
  int        nt = 0;
  int        status;
   
  // watch out!
  // a new sequence starts with an integer,
  // the chainname can also be an interger, and it appears at the begining of line
  // a newsequence line starts with an integer and it has at least 4 tokens
  // a chainname line has only two tokens
  //
  esl_sprintf(&tmp, line);
  if ((status = esl_strtok(&tmp, " \t", &tok)) == eslOK) {
    if (esl_str_IsInteger(tok)) newseq = TRUE;
    seqid = atoi(tok);
  }

  if (newseq) {
    nt ++;

    // second token is the type of sequence (dna, rna, peptide...)
    if ( (status = esl_strtok(&tmp, " \t", &tok)) == eslOK) {
      esl_sprintf(&seqtype, tok); nt ++;
    }

    // if there is a third token, this is definitively a new sequence
    if ( (status = esl_strtok(&tmp, " \t", &tok)) == eslOK) {
      if (chain->seq)  free(chain->seq);  chain->seq  = NULL;
      if (chain->name) free(chain->name); chain->name = NULL;
      
      pdbx->nsq ++;
      chain->seqid = seqid;
      esl_sprintf(&chain->seqtype, seqtype);
      if (verbose) printf("new sequence ");
  
      nt ++;
    }
    else return FALSE;
  
    // if there is a fifth token, that means the sequence fits in one tmp
    while ( (status = esl_strtok(&tmp, " \t", &tok)) == eslOK) {
      nt ++;

      // nt=3,4 not of interest here.
      // nt=5   seq (we are not interested on it)

      if (nt == 5) {
	isshort = TRUE;
	if (verbose) printf("short? %d\n", isshort);
      }
      
      if (nt == 6) { // the seqcan fitted here too
	esl_strdup(tok, -1, &chain->seq);
	if (verbose) printf("seq %s\n", chain->seq);
      }
      if (nt == 7) { // the chain fitted here too
	// but we are not adding the chainname from here
	// we are passing it back so we know the canonical sequence is all in,
	// but the chain name will not be stored.
	esl_strdup(tok, -1, &chain->name);
     }
    }
    if (nt < 5 && verbose) printf("short? %d\n", isshort);
  }

  *ret_isshort = (newseq)? isshort : *ret_isshort;
  
  if (seqtype) free(seqtype);
  return newseq;
}

static void
chunk_is_sequence_short(char *line, char **ret_seq, char **ret_chainname, int verbose)
{
  char *tok;
  int   status;

  // if chainname exists, we are done
  if (*ret_chainname) return;

  // else if seqcan exists, first token is chain,
  // otherwise first token is the consensus sequence
  //
  if ( (status = esl_strtok(&line, " \t", &tok)) == eslOK) {
    if      (*ret_seq == NULL) esl_strdup(tok, -1, ret_seq);
    else if (ret_chainname)    esl_strdup(tok, -1, ret_chainname);
  }
  if (*ret_seq && verbose) printf("seq %s\n", *ret_seq);

}

static int
chunk_is_sequence_start(char *line, char **ret_seq, SEQSTATUS *ret_seqstatus)
{
  SEQSTATUS  seqstatus = *ret_seqstatus;
  char      *tok;
  int        isstart = FALSE;
  int        status;

  if ( (status = esl_strtok(&line, " \t", &tok)) == eslOK) {
    
    // a new sequence
    if (strlen(tok) > 1 && tok[0] == ';') {

      isstart = TRUE;
      // _entity_poly.pdbx_seq_one_letter_code
      // we don't save this sequence
      if (seqstatus == SEQ_NONE) seqstatus = SEQ_START;
      
      // _entity_poly.pdbx_seq_one_letter_code_can
      // we save the canonical sequence
      else if (seqstatus == SEQ_END)  {
	seqstatus = SEQCAN_START;
	
	// ignore the starting ';'
	esl_strdup(tok+1, -1, ret_seq);
      }
      
      // otherwise, you should not be here
      else  esl_fatal("fatal at sequence_start(): seqstatus %d", seqstatus);
    }
    
    // the weird case with lots of (UNK)
    // the sequence fits all here and we end
    else if (strlen(tok) > 1 && tok[0] != ';' && seqstatus == SEQ_END) {
      esl_strdup(tok, -1, ret_seq);
      seqstatus = SEQCAN_END;
    }
    
  }
  
  *ret_seqstatus = seqstatus;
  return isstart;
}

static int
chunk_is_sequence_cont(char *line, char **ret_seq, SEQSTATUS *ret_seqstatus)
{
  SEQSTATUS  seqstatus = *ret_seqstatus;
  char      *tok;
  int        iscont = FALSE;
  int        status;

  if (seqstatus == SEQCAN_START || seqstatus == SEQCAN_CONT) {
    if ( (status = esl_strtok(&line, " \t\n", &tok)) == eslOK) {
      // a cont sequence
      if (tok[0] != ';') {
	
	iscont = TRUE;
	
	// SEQ:     _entity_poly.pdbx_seq_one_letter_code 
	// SEQCANL: _entity_poly.pdbx_seq_one_letter_code_can (the one we want)
	//
	esl_strcat(ret_seq, -1, tok, -1);
	seqstatus = SEQCAN_CONT;
      }
    }
  }
  
  *ret_seqstatus = seqstatus;
  return iscont;
}

static int
chunk_is_sequence_end(char *line, SEQSTATUS *ret_seqstatus)
{
  SEQSTATUS   seqstatus = *ret_seqstatus;
  char       *tok;
  int         isend = FALSE;
  int         status;

  if (seqstatus == SEQ_NONE) esl_fatal("fatal at sequence_end()");
  
  // a lone ';' character
  if ( (status = esl_strtok(&line, " \t", &tok)) == eslOK) {
    if (strlen(tok) == 1 && tok[0] == ';') {
      isend = TRUE;
      if (seqstatus == SEQ_START    || seqstatus == SEQ_CONT)    seqstatus = SEQ_END;
      if (seqstatus == SEQCAN_START || seqstatus == SEQCAN_CONT) seqstatus = SEQCAN_END;  
    }    
 
  }
  *ret_seqstatus = seqstatus;
  return isend; 
}

static int
chunk_is_chainname(char *line, char **ret_name)
{
  char *tok;
  int   isname = FALSE;
  int   status;
  
  // first token is chain
  if ( (status = esl_strtok(&line, " \t", &tok)) == eslOK) {
    if (tok[0] != ';') {
      if (*ret_name) esl_strdup(tok, -1, ret_name);
      isname = TRUE;
    }
  }

   return isname;
}


// this loop gives us the correspondence between
// a chain (_struct_asym.id) and a sequence (_struct_asym.entity_id)
//
// # 
// loop_
// _struct_asym.id 
// _struct_asym.pdbx_blank_PDB_chainid_flag 
// _struct_asym.pdbx_modified 
// _struct_asym.entity_id 
// _struct_asym.details 
// A N N 1 ? 
// B N N 2 ? 
// C N N 2 ? 
// D N N 2 ? 
//
static int
parse_sasymfield(char *p, esl_pos_t n, int nsl, char **sasymlabel, int *ret_nsaf, int **ret_usesasymf, int verbose)
{
  int  nsaf      = *ret_nsaf;
  int *usesasymf = *ret_usesasymf;
  int  s;
  int  status;

  if (nsaf == 0) ESL_ALLOC  (usesasymf, sizeof(int));
  else           ESL_REALLOC(usesasymf, sizeof(int) * (nsaf+1));
  
  usesasymf[nsaf] = -1; // sasym field not used by default
  for (s = 0; s < nsl; s ++) {
    if  (esl_memstrpfx(p, n, sasymlabel[s])) {
      usesasymf[nsaf] = s;
      if (verbose) printf ("%s is struct_saym field %d\n", sasymlabel[s], nsaf);
    }
  }
  
  *ret_nsaf      = ++nsaf;
  *ret_usesasymf = usesasymf;
  return eslOK;

 ERROR:
  if (usesasymf) free(usesasymf);
  return status;
}

static int
parse_sasym(PDBX *pdbx, char *p, esl_pos_t n, int nsf, int *usesasymfield, int verbose)
{
  struct chain_s  *chain = NULL;
  char           **field = NULL;
  esl_pos_t        i;
  esl_pos_t        len     = 0;
  esl_pos_t        salloc  = n + 1;
  int              idx     = 0;
  int              which   = -1;
  int              isblank = FALSE;
  int              f;
  int              c;
  int              status;

  chain = rview_chain_Create();
  ESL_ALLOC(field,    sizeof(char *) * nsf);
  ESL_ALLOC(field[0], sizeof(char  ) * nsf * salloc);
  for (f = 1; f < nsf; f ++) field[f] = field[f-1] + salloc;
  
  for (i = 0; i < n; i++) {
    if (! isspace(p[i])) {
      isblank = FALSE;
      field[idx][len++] = p[i];
    }
    else {
      if (isblank == FALSE) { field[idx][len] = '\0';  idx ++; }
      isblank = TRUE;
      len     = 0;
    }
  }
  if (idx != nsf) esl_fatal("wrong number of struct_asym fields %d should be %d\n", idx, nsf);  

  // reverse the hard-coded labels into
  // chain variables (name and seqid)
  //
  for (f = 0; f < nsf; f++) {
    if (usesasymfield[f] ==  0) esl_sprintf(&chain->name, field[f]);
    if (usesasymfield[f] ==  1) chain->seqid = atoi(field[f]);
  }
  
  // Look at pdbx,
  // if the last chain associated to this sequence has no name (chain->name == NULL),
  // just add the name. Otherwise add a new chain.
  //
  for (c = 0; c < pdbx->nch; c ++)
    if (chain->seqid == pdbx->chain[c].seqid) which = c;
  if (which == pdbx->nch) esl_fatal("you should not be here. No sequence has been found with id %d\n", chain->seqid);
 
 // If which < 0, the chain is not associated to a sequence.
  // It could be a co-factor, or others
  if (which >=  0) {

    // there is only one chain for this sequence, just add the name
    if (pdbx->chain[which].name == NULL && pdbx->chain[which].seq) {
      esl_sprintf(&(pdbx->chain[which].name), chain->name);
      if (verbose)  rview_chain_Write(stdout, &pdbx->chain[which], verbose);
    }
    // add one more chain with the same sequence, and the new name
    else {
      esl_sprintf(&chain->seq,     pdbx->chain[which].seq);
      esl_sprintf(&chain->seqtype, pdbx->chain[which].seqtype);
      rview_pdbx_AddChain(pdbx, chain, which+1, verbose);
     }
  }

  rview_chain_Destroy(1, chain);
  free(field[0]);
  free(field);
  return eslOK;

 ERROR:
  if (field[0]) free(field[0]);
  if (field)    free(field);
  if (chain)    rview_chain_Destroy(1, chain);
  return status;
}

// in ATOMS amino acids are given in the three letter alphabet.
// this function returns its one-letter traslation.
static int
res_oneletter(char **ret_x, char *res)
{
  char *x = NULL;

  if      (strlen(res) == 1)                  esl_sprintf(&x, res);
  else if (strlen(res) == 3) {    
    if      (esl_strcmp(res, "ALA") == eslOK) esl_sprintf(&x, "A");
    else if (esl_strcmp(res, "ARG") == eslOK) esl_sprintf(&x, "R");
    else if (esl_strcmp(res, "ASN") == eslOK) esl_sprintf(&x, "N");
    else if (esl_strcmp(res, "ASP") == eslOK) esl_sprintf(&x, "D");
    else if (esl_strcmp(res, "ASX") == eslOK) esl_sprintf(&x, "B");
    else if (esl_strcmp(res, "CYS") == eslOK) esl_sprintf(&x, "C");
    else if (esl_strcmp(res, "GLU") == eslOK) esl_sprintf(&x, "E");
    else if (esl_strcmp(res, "GLN") == eslOK) esl_sprintf(&x, "Q");
    else if (esl_strcmp(res, "GLX") == eslOK) esl_sprintf(&x, "Z");
    else if (esl_strcmp(res, "GLY") == eslOK) esl_sprintf(&x, "G");
    else if (esl_strcmp(res, "HIS") == eslOK) esl_sprintf(&x, "H");
    else if (esl_strcmp(res, "ILE") == eslOK) esl_sprintf(&x, "I");
    else if (esl_strcmp(res, "LEU") == eslOK) esl_sprintf(&x, "L");
    else if (esl_strcmp(res, "LYS") == eslOK) esl_sprintf(&x, "K");
    else if (esl_strcmp(res, "MET") == eslOK) esl_sprintf(&x, "M");
    else if (esl_strcmp(res, "PHE") == eslOK) esl_sprintf(&x, "F");
    else if (esl_strcmp(res, "PRO") == eslOK) esl_sprintf(&x, "P");
    else if (esl_strcmp(res, "SER") == eslOK) esl_sprintf(&x, "S");
    else if (esl_strcmp(res, "THR") == eslOK) esl_sprintf(&x, "T");
    else if (esl_strcmp(res, "TRP") == eslOK) esl_sprintf(&x, "W");
    else if (esl_strcmp(res, "TYR") == eslOK) esl_sprintf(&x, "Y");
    else if (esl_strcmp(res, "VAL") == eslOK) esl_sprintf(&x, "V");
    else if (esl_strcmp(res, "UNK") == eslOK) esl_sprintf(&x, "X");
    else                                      esl_fatal("cannot find amino acid %s\n", res);
  }
  else                                        esl_fatal("cannot identify this residue %s\n", res);
  
  *ret_x = x;
  return eslOK;
}

static void
writep(char *p, esl_pos_t n)
{ int i;
  for (i = 0; i < n; i ++) printf("%c", p[i]);
  printf("\n");
}
