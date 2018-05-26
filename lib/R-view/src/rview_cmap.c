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
  PDBX  *pdbx  = NULL;
  CLIST *clist = NULL;
  int    nct   = 0;
  int    status;

  status = rview_ReadPDBxfile(pdbxfile, &pdbx, errbuf, verbose);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "rview_ReadPDBxfile() failed");
  
  status = rview_CreateContacts(pdbx, &nct, &clist, errbuf, verbose);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "rview_CreateContacts() failed");
  
  if (pdbx) rview_pdbx_Destroy(pdbx);
  
  *ret_nct   = nct;
  *ret_clist = clist; 
  return eslOK;

 ERROR:
  if (pdbx) rview_pdbx_Destroy(pdbx);
  if (clist) CMAP_FreeCList(clist);
  return status;
}

