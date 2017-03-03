/*  homolog - given two seqs, calculate homology
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
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sqio.h"
#include "esl_tree.h"

#include "hmmer.h"

#include "homology.h"


int
homol_RunPHMMER(int n, int m, const ESL_MSA *msa, float *ret_sc, float *ret_eval, char *errbuf, int verbose)
{
  char       phmmertbl[16] = "esltmpXXXXXX"; 
  char       phmmerout[16] = "esltmpXXXXXX"; 
  char       sqfile[16] = "esltmpXXXXXX";
  char       dbfile[16] = "esltmpXXXXXX"; 
  FILE      *phmmerfp = NULL;
  FILE      *sqfp = NULL;
  FILE      *dbfp = NULL;
  char      *s = NULL;
  char      *args = NULL;
  ESL_SQ    *sq1 = NULL;
  ESL_SQ    *sq2 = NULL;
  float      sc;
  float      eval;
  int        status;
  
  if ("EVOHMMDIR" == NULL)               { status = eslENOTFOUND; goto ERROR; }
  if ((s = getenv("EVOHMMDIR")) == NULL) { status = eslENOTFOUND; goto ERROR; }

  /* phmmer tblout file */
  if ((status = esl_tmpfile_named(phmmertbl, &phmmerfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create phmmertbl");
  fclose(phmmerfp);
  if ((status = esl_tmpfile_named(phmmerout, &phmmerfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create phmmerout");
  fclose(phmmerfp);

  /* put seqs in fasta format */
  if ((status = esl_tmpfile_named(sqfile, &sqfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create sqfile");
  esl_sq_FetchFromMSA(msa, n, &sq1);
  esl_sqio_Write(sqfp, sq1, eslSQFILE_FASTA, FALSE);
  fclose(sqfp);
  if ((status = esl_tmpfile_named(dbfile, &dbfp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create dbfile");
  esl_sq_FetchFromMSA(msa, m, &sq2);
  esl_sqio_Write(dbfp, sq2, eslSQFILE_FASTA, FALSE);
  fclose(dbfp);
 
  /* run phmmer */
  esl_sprintf(&args, "%s/lib/hmmer4/src/programs/phmmer --tbl %s -o %s %s %s\n", s, phmmertbl, phmmerout, sqfile, dbfile);
  system(args);
  if (args) free(args); args = NULL;
 
  /* parse phmmer tblout */
  if ((status = homol_ParsePHMMERtblout(phmmertbl, &sc, &eval, errbuf, verbose)) != eslOK) ESL_XFAIL(status, errbuf, "failed to parse phmmer tblout");

  *ret_sc   = sc;
  *ret_eval = eval;

  remove(phmmertbl);
  remove(phmmerout);
  remove(sqfile);
  remove(dbfile);
  esl_sq_Destroy(sq1);
  esl_sq_Destroy(sq2);
  return eslOK;

 ERROR:
  if (args) free(args);
  if (sq1)  esl_sq_Destroy(sq1);
  if (sq2)  esl_sq_Destroy(sq2);
  remove(phmmertbl);
  remove(phmmerout);
  remove(sqfile);
  remove(dbfile);
  return status;
}

int
homol_ParsePHMMERtblout(char *tblout, float *ret_sc, float *ret_eval, char *errbuf, int verbose)
{
  FILE  *fp = NULL;
  char  *buf = NULL;           /* growable buffer for esl_fgets()      */
  char  *s;		
  char  *sqname1;
  char  *sqname2;
  char  *sqacc1;
  char  *sqacc2;
  char  *ceval;
  char  *cscore;
  float  sc   = -eslINFINITY;
  float  eval = eslINFINITY;
  int    n = 0;
  int    status;
  
  if ((fp = fopen(tblout, "r"))  == NULL) ESL_XFAIL(eslFAIL, errbuf, "failed to open %s\n", tblout); 

  while (esl_fgets(&buf, &n, fp) == eslOK)
    {
      if (*buf == '#') continue; /* skip comments */
      
      s = buf;                      
      esl_strtok(&s, " \t", &sqname1);     /* sq1 name   = 1st token online */
      esl_strtok(&s, " \t", &sqacc1);      /* sq1 acc    = 2st token online */
      esl_strtok(&s, " \t", &sqname2);     /* sq2 name   = 3st token online */
      esl_strtok(&s, " \t", &sqacc2);      /* sq2 acc    = 4st token online */
      esl_strtok(&s, " \t", &ceval);       /* eval       = 5st token online */
      esl_strtok(&s, " \t", &cscore);      /* score      = 6st token online */

      if (!esl_str_IsReal(cscore) || !esl_str_IsReal(ceval)) ESL_XFAIL(eslFAIL, errbuf, "failed to parse eval %s or score %s\n", ceval, cscore);
      sc   = atof(cscore);
      eval = atof(ceval);
     }  
  
  *ret_sc   = sc;
  *ret_eval = eval;

  fclose(fp);
  free(buf);

  return eslOK; 
  
 ERROR:
  if (buf) free(buf);
  return status;
}



/*****************************************************************
 * @LICENSE@
 *****************************************************************/
