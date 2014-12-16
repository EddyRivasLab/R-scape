/*  mya - get distance in millions of years between two species
 * 
 * Contents:
 *   1. Miscellaneous functions for msatree
 *   2. Unit tests
 *   3. Test driver
 *   4. License and copyright 
 *
 * ER, Fri Oct 21 10:26:52 EDT 2011 [Janelia] 
 * SVN $Id:$
 */

#include "p7_config.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <math.h>
#include <float.h>
	
#include "easel.h"

#include "hmmer.h"

#include "mya.h"
#include "fetchpfamdb.h"
#include "ssifile.h"


int
mya_BetweenSpeciesFromFile(char *distfile, char *tx1, char *tx2, float *ret_mya, int verbose)
{
  char   *ssidistfile = NULL;
  char  **taxalist = NULL;
  FILE   *fp = NULL;
  float   mya = 0.0;
  int     ntaxa;
  int     x;
  int     status;

  if (distfile == NULL)  { *ret_mya = -1.0; return eslOK; }
  
  /* sometimes ncbi cannot assign a species id to a given Pfam seqid.
   */
  if (tx1 == NULL || tx2 == NULL) { *ret_mya = -1.0; return eslOK; }
  
 /* Create an ssi index of spfile.
   * spfile contains the info to associate species to Pfam sq names */
  status = ssi_distfile_Create(distfile, &ssidistfile, &taxalist, &ntaxa);
  if (status != eslOK) { printf("ssi_dist_Create failed for file %s\n", distfile); status = eslFAIL; goto ERROR; }
  if (verbose) printf("ssifile %s created\n", ssidistfile);

  if ((fp = fopen(distfile, "r")) == NULL) { printf("failed to open %s\n", distfile); status = eslFAIL; goto ERROR;  }

  ssi_distfile_GetDistance(ssidistfile, taxalist, ntaxa, tx1, tx2, &mya);  
  if (verbose) printf(" D(%s, %s) = %f\n", tx1, tx2, mya);
   
  fclose(fp);
  remove(ssidistfile);
  for (x = 0; x < ntaxa; x ++) free(taxalist[x]);
  free(taxalist);

  *ret_mya = mya;
  return eslOK;

 ERROR:
  for (x = 0; x < ntaxa; x ++) if (taxalist[x]) free(taxalist[x]);
  if (taxalist) free(taxalist);
  return status;
}

int
mya_BetweenSpeciesFromSSIFile(char *ssidistfile, char **taxalist, int ntaxa, char *tx1, char *tx2, float *ret_mya, int verbose)
{
  float  mya = 0.0;
  int    status;

  if (ssidistfile == NULL)  { *ret_mya = -1.0; return eslOK; }
  
  /* sometimes ncbi cannot assign a species id to a given Pfam seqid.
   */
  if (tx1 == NULL || tx2 == NULL) { *ret_mya = -1.0; return eslOK; }
  
  status = ssi_distfile_GetDistance(ssidistfile, taxalist, ntaxa, tx1, tx2, &mya);
  if (status != eslOK) goto ERROR;
  if (verbose) printf(" D(%s, %s) = %f\n", tx1, tx2, mya);
   
  *ret_mya = mya;
  return eslOK;

 ERROR:
  return status;
}

int
mya_BetweenSpecies(char *tx1, char *tx2, float *ret_mya, char *errbuf, int verbose)
{
  char            distfile[16] = "esltmpXXXXXX"; /* temp distfile */
  char           *args = NULL;
  char           *s;
  char           *tok;
  char           *tok1, *tok2;
  FILE           *distfp = NULL;
  ESL_FILEPARSER *efp = NULL;
  float           mya = -1.0;
  int             status;
  
  /* sometimes ncbi cannot assign a species id to a given Pfam seqid.
   */
  if (tx1 == NULL || tx2 == NULL) { *ret_mya = -1.0; return eslOK; }
  
  if ("EVOHMMDIR" == NULL)               { status = eslENOTFOUND; goto ERROR; }
  if ((s = getenv("EVOHMMDIR")) == NULL) { status = eslENOTFOUND; goto ERROR; }
  
  /* Special case: comparison to itself. 
   * timetree cannot deal with returning a 0 distance 
   */
  if (esl_strcmp(tx1, tx2) == 0) { *ret_mya = 0.0; return eslOK; }

  if ((status = esl_tmpfile_named(distfile, &distfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create distfile");
  fclose(distfp);

  /* run script */
  esl_sprintf(&args, "%s/scripts/fetchTimeTreeData.pl -ta '%s' -tb '%s' > %s\n", s, tx1, tx2, distfile);
  system(args);
  
 /* get mya from distfile */
  if (esl_fileparser_Open(distfile, NULL, &efp) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to open file %s", distfile);

   while (esl_fileparser_NextLine(efp) == eslOK)
     {
       if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse file %s", distfile);

       if (strcmp(tok, "Mean") == 0) {
	 if (esl_fileparser_GetTokenOnLine(efp, &tok1, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse file %s", distfile);  
	 if (esl_fileparser_GetTokenOnLine(efp, &tok2, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "failed to parse file %s", distfile);  
	 if (strcmp(tok1, "is") == 0) {
	   mya = atof(tok2);
	 }
       }
     }
   esl_fileparser_Close(efp);

   if (verbose) printf(" MYA(%s,%s) = %f\n", tx1, tx2, mya);
  
  *ret_mya = mya;  

  remove(distfile);
  free(args);
  return eslOK;

 ERROR:
  if (args != NULL) free(args);
  remove(distfile);
  return status;
}

int
mya_BetweenLevels(int n, int m, SPECIES *SPE, float *ret_mya, char *errbuf, int verbose)
{
  float  mya = -1.0;
  int    t1, t2;
  int    status;
 
  /* first try using the especies names */
  status = mya_BetweenSpecies(SPE->spname[n], SPE->spname[m], &mya, errbuf, verbose); if (status != eslOK) goto ERROR;
  if (mya >= 0) { *ret_mya = mya; return eslOK; }
 
  t1 = NTAXALEVELS-1;
  t2 = NTAXALEVELS-1;
  status = mya_BetweenSpecies(SPE->spname[n],     SPE->parent[m][t2], &mya, errbuf, verbose); if (status != eslOK) goto ERROR; if (mya >= 0) { *ret_mya = mya; return eslOK; }
  status = mya_BetweenSpecies(SPE->parent[n][t1], SPE->spname[m],     &mya, errbuf, verbose); if (status != eslOK) goto ERROR; if (mya >= 0) { *ret_mya = mya; return eslOK; }
  
  for (t1 = NTAXALEVELS-1; t1 >= 0; t1 --) 
    for (t2 = NTAXALEVELS-1; t2 >= t1; t2 --) {
      if (esl_strcmp(SPE->parent[n][t1], SPE->parent[m][t2]) == 0) { *ret_mya = mya; return eslOK; }
      
      /* if you have not reached the lca node yet */
      status = mya_BetweenSpecies(SPE->parent[n][t1], SPE->parent[m][t2], &mya, errbuf, verbose); if (status != eslOK) goto ERROR; 
      if (mya >= 0) { *ret_mya = mya; return eslOK; } 
    }
 
  *ret_mya = mya;  

  return eslOK;

 ERROR:
  return eslFAIL;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
