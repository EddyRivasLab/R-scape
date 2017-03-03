/*  muscle -- function to integrate MUSCLE msa output
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
MUSCLE_Align(const ESL_MSA *msa, ESL_MSA **ret_musclemsa, char *errbuf, int verbose)
{
  ESL_MSAFILE *afp = NULL;
  FILE         *fp = NULL;
  ESL_MSA      *musclemsa = NULL;
  char          tmpinfile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char         *args = NULL;
  char         *s = NULL;
  int           status;
  
  /* MUSCLE input and output in AFA format */
  if ((status = esl_tmpfile_named(tmpinfile,  &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create input file");
  if ((status = esl_msafile_Write(fp, (ESL_MSA *)msa, eslMSAFILE_AFA)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write AFA file\n");
  fclose(fp);
  
  /* run MUSCLE */
  if ("MUSCLEDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("MUSCLEDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpoutfile, &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  //esl_sprintf(&args, "%s/muscle -in %s -out %s", s, tmpinfile, tmpoutfile);
  esl_sprintf(&args, "%s/muscle3.8.31_i86linux64 -in %s -out %s", s, tmpinfile, tmpoutfile);
  system(args);
  fclose(fp);
  
  /* convert output to msa */
  if (esl_msafile_Open(NULL, tmpoutfile, NULL, eslMSAFILE_AFA, NULL, &afp) != eslOK) ESL_XFAIL(status, errbuf, "Failed to open AFA file\n");
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);
  afp->format = eslMSAFILE_AFA;
  if (esl_msafile_Read(afp, &musclemsa) != eslOK) ESL_XFAIL(status, errbuf, "Failed to read AFA file\n");
  if (status != eslOK) esl_msafile_ReadFailure(afp, status);
  esl_msafile_Close(afp);
  
  remove(tmpinfile);
  remove(tmpoutfile);
  
  *ret_musclemsa = musclemsa;

  if (args != NULL) free(args);
  return eslOK;
  
 ERROR:
  remove(tmpinfile);
  remove(tmpoutfile);
  if (args != NULL) free(args);
  return status;  

}
