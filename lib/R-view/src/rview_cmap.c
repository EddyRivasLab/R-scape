/* rview_cmap - the R-view IPA
 *
 */
#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_fileparser.h"

#include "rview_contacts.h"
#include "rview_cmap.h"
#include "rview_pdbfile.h"


int rview_ContactMap(char *pdbxfile, int *ret_nct, CLIST **ret_clist, double maxD, int minL, char *errbuf, int verbose)
{
  struct chain_s *chain = NULL;
  CLIST          *clist = NULL;
  int              nchain = 0;
  int              nct    = 0;
  int              status;

  status = rview_ReadPDBxfile(pdbxfile, &nchain, &chain, errbuf, verbose);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "rview_ReadPDBxfile() failed");
  
  status = rview_CreateContacts(nchain, chain, &nct, &clist, errbuf, verbose);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "rview_CreateContacts() failed");
  
  if (chain) rview_chain_Destroy(chain);
  
  *ret_nct   = nct;
  *ret_clist = clist; 
  return eslOK;

 ERROR:
  if (chain) rview_chain_Destroy(chain);
  if (clist) CMAP_FreeCList(clist);
  return status;
}

