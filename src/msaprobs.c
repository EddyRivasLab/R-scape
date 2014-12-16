/*  msaprobs -- function to integrate MSAProbs msa output
 *
 * ER, Mon Feb 24 15:44:55 EST 2014 [Janelia] 
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
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "muscle.h"

int
MSAProbs_Align(const ESL_MSA *msa, ESL_MSA **ret_msaprobsmsa, char *errbuf, int verbose)
{
  ESLX_MSAFILE *afp = NULL;
  FILE         *fp = NULL;
  ESL_MSA      *msaprobsmsa = NULL;
  char          tmpinfile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char         *args = NULL;
  char         *s = NULL;
  int           status;
  
  /* MSAProbs input and output in AFA format */
  if ((status = esl_tmpfile_named(tmpinfile,  &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create input file");
  if ((status = eslx_msafile_Write(fp, (ESL_MSA *)msa, eslMSAFILE_AFA)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write AFA file\n");
  fclose(fp);
  
  /* run MSAProbs */
  if ("MSAPROBSDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("MSAPROBSDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpoutfile, &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  esl_sprintf(&args, "%s/MSAProbs/msaprobs %s > %s", s, tmpinfile, tmpoutfile);
  system(args);
  fclose(fp);
  
  /* convert output to msa */
  if (eslx_msafile_Open(NULL, tmpoutfile, NULL, eslMSAFILE_AFA, NULL, &afp) != eslOK) eslx_msafile_OpenFailure(afp, status);
  afp->format = eslMSAFILE_AFA;
  if (eslx_msafile_Read(afp, &msaprobsmsa) != eslOK) eslx_msafile_ReadFailure(afp, status);
  eslx_msafile_Close(afp);
  
  remove(tmpinfile);
  remove(tmpoutfile);
  
  *ret_msaprobsmsa = msaprobsmsa;

  if (args != NULL) free(args);
  return eslOK;
  
 ERROR:
  remove(tmpinfile);
  remove(tmpoutfile);
  if (args != NULL) free(args);
  return status;  

}
