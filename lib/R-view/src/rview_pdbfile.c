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
#include "esl_fileparser.h"
#include "esl_sq.h"

#include "rview_pdbfile.h"

int
rview_ReadPDBxfile(char *pdbxfile, int *ret_nchain, struct chain_s **ret_chain, char *errbuf, int verbose)
{
  ESL_FILEPARSER  *fp    = NULL;
  struct chain_s  *chain = NULL;
  char            *tok;
  int              nchain = 0;
  int              status;

  // Read the pdbfile
  //
  // For each chain,
  // we are going to extract the sequence (chain->seq) and the information about all the residues and their atoms (chain->res)
  //
  if (esl_fileparser_Open(pdbxfile, NULL, &fp) != eslOK)  ESL_XFAIL(eslFAIL, errbuf, "file open failed");
  esl_fileparser_SetCommentChar(fp, '#');
  
  while (esl_fileparser_NextLine(fp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(fp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse token from file %s", pdbxfile);
      
    }
  esl_fileparser_Close(fp);

  // create the map between chain->seq and the sequence described by the atoms

  *ret_nchain = nchain;
  return eslOK;

 ERROR:
  if (chain) rview_chain_Destroy(chain);
  return status;
}

void
rview_atom_Init(ATOM *atom)
{
  atom->chain = NULL;
  atom->idx   = -1;
  atom->type  = ATYPE_NONE;
}

void
rview_atom_Destroy(ATOM *atom)
{
  if (atom->chain) free(atom->chain);
}

void
rview_chain_Init(struct chain_s *chain)
{
  chain->L        = 0;
  chain->seq      = NULL;
  
  chain->nr       = 0;
  chain->res      = NULL;
  
  chain->atom_map = NULL; 
}

int
rview_chain_AddRes(struct chain_s *chain)
{
  int status;
  if (chain->nr == 0) ESL_ALLOC  (chain->res, sizeof(RES));
  else                ESL_REALLOC(chain->res, sizeof(RES) * (chain->nr+1));

  chain->nr ++;
  return eslOK;

 ERROR:
  return eslFAIL;
}

void
rview_chain_Destroy(struct chain_s *chain)
{
  int n;

  if (chain->seq)    esl_sq_Destroy(chain->seq);
  if (chain->resseq) esl_sq_Destroy(chain->resseq);
  if (chain->res) {
    for (n = 0; n < chain->nr; n ++) rview_res_Destroy(&chain->res[n]);
    free(chain->res);
  }
  if (chain->atom_map) free(chain->atom_map);
}

int
rview_res_AddAtom(RES *res)
{
  int status;
  if (res->na == 0) ESL_ALLOC  (res->atom, sizeof(ATOM));
  else              ESL_REALLOC(res->atom, sizeof(ATOM) * (res->na+1));

  res->na ++;
  return eslOK;
 ERROR:
  return eslFAIL;
}

void
rview_res_Init(RES *res)
{
  res->chain = NULL;

  res->from_atomidx = -1;
  res->resnum       = -1;

  res->na    = 0;
  res->atom  = NULL;
}

void
rview_res_Destroy(RES *res)
{
  int n;
  
  if (res->chain) free(res->chain);
  if (res->atom) {
    for (n = 0; n < res->na; n ++) rview_atom_Destroy(&res->atom[n]);
    free(res->atom);
  }
}

