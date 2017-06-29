#include "esl_msa.h"
/* fastSP.c -- function to integrate FastSP output
 *
 * ER, Wed Feb 26 10:15:35 EST 2014 [Janelia] 
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
#include "esl_fileparser.h"
#include "esl_msa.h"
#include "esl_msafile.h"

#include "hmmer.h"

#include "fastsp.h"

static double calculateF(double sen, double ppv);
static int    fastsp_parse(char *fastspfile, int *ret_tph, int *ret_trueh, int *ret_foundh, int *ret_cac, int *ret_ac, double *ret_sen, double *ret_ppv, 
			   double *ret_F, double *ret_TC, int flip, char *errbuf, int verbose);

int
FastSP_Run(const ESL_MSA *msar, const ESL_MSA *msae, int *ret_tph, int *ret_trueh, int *ret_foundh, int *ret_cac, int *ret_ac,
	   double *ret_sen, double *ret_ppv, double *ret_F, double *ret_TC, int lcu_r, int lcu_e, char *errbuf, int verbose)
{
  FILE         *fp = NULL;
  char          fastspv[17]    = "FastSP_1.6.0.jar";
  char          tmpmsar[16]    = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpmsae[16]    = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char         *args = NULL;
  char         *s = NULL;
  int           trueh  = 0;
  int           foundh = 0;
  int           tph    = 0;
  int           cac    = 0;     // correctly aligned columns
  int           ac     = 0;     // aligned columns in reference alignment
  double        sen    = 0.0;
  double        ppv    = 0.0;
  double        F      = 0.0;
  double        TC     = 0.0;
  int           status;
  

  /* FastSP reference and estimated MSAs in AFA format */
  if ((status = esl_tmpfile_named(tmpmsar,  &fp))                        != eslOK) ESL_XFAIL(status, errbuf, "failed to create msar file");
  if ((status = esl_msafile_Write(fp, (ESL_MSA *)msar, eslMSAFILE_AFA)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write AFA file\n");
  fclose(fp);
  if (msae) {
    if ((status = esl_tmpfile_named(tmpmsae,  &fp))                        != eslOK) ESL_XFAIL(status, errbuf, "failed to create msae file");
    if ((status = esl_msafile_Write(fp, (ESL_MSA *)msae, eslMSAFILE_AFA)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write AFA file\n");
    fclose(fp);
  }
  
  /* run FastSP */
  if ("FASTSPDIR" == NULL)               ESL_XFAIL(eslENOTFOUND, errbuf, "FASTSPDIR not found");
  if ((s = getenv("FASTSPDIR")) == NULL) ESL_XFAIL(eslENOTFOUND, errbuf, "FASTSPDIR not found");
  if ((status = esl_tmpfile_named(tmpoutfile, &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  // option -ml  considers lower case in the estimated msa as unaligned
  // option -mlr considers lower case in the reference msa as unaligned (-ml and -mlr can be used together)
  // option -c .- tells the program that both . and - should be considered gaps

  if (lcu_e && lcu_r) 
    esl_sprintf(&args, "java -jar -Xmx2048m %s/%s -r  %s -e %s -c .-~_ -ml  -mlr &> %s", s, fastspv, tmpmsar, (msae)? tmpmsae:tmpmsar, tmpoutfile); // direct stderr to file too
  else if (lcu_e)
    esl_sprintf(&args, "java -jar -Xmx2048m %s/%s -r  %s -e %s -c .-~_ -ml       &> %s", s, fastspv, tmpmsar, (msae)? tmpmsae:tmpmsar, tmpoutfile); // direct stderr to file too
  else if (lcu_r)
    esl_sprintf(&args, "java -jar -Xmx2048m %s/%s -r  %s -e %s -c .-~_      -mlr &> %s", s, fastspv, tmpmsar, (msae)? tmpmsae:tmpmsar, tmpoutfile); // direct stderr to file too
  else 
    esl_sprintf(&args, "java -jar -Xmx2048m %s/%s -r  %s -e %s -c .-~_           &> %s", s, fastspv, tmpmsar, (msae)? tmpmsae:tmpmsar, tmpoutfile); // direct stderr to file too
  system(args);
 if ((status = fastsp_parse(tmpoutfile, &tph, &trueh, &foundh, &cac, &ac, &sen, &ppv, &F, &TC, FALSE, errbuf, verbose)) != eslOK) 
    ESL_XFAIL(status, errbuf, "failed to extract F");
  fclose(fp);
   
  if (msae == NULL) {
    tph    = 0;
    foundh = 0;
    cac    = 0;
    ac     = 0;
    sen    = 0.;
    ppv    = 0.;
    F      = 0.;
    TC     = 0.;
  }

  if (ret_tph)    *ret_tph    = tph;
  if (ret_trueh)  *ret_trueh  = trueh;
  if (ret_foundh) *ret_foundh = foundh;
  if (ret_cac)    *ret_cac    = cac;
  if (ret_ac)     *ret_ac     = ac;
  if (ret_sen)    *ret_sen    = sen;
  if (ret_ppv)    *ret_ppv    = ppv;
  if (ret_F)      *ret_F      = F;
  if (ret_TC)     *ret_TC     = TC;
 
  remove(tmpmsar);
  remove(tmpmsae);
  remove(tmpoutfile);
   
  if (args != NULL) free(args);
  return eslOK;
  
 ERROR:
  remove(tmpmsar);
  remove(tmpmsae);
  remove(tmpoutfile);
  if (args != NULL) free(args);
  return status;  

}

int
FastSP_Benchmark(FILE *benchfp, char *msaname, char *method, ESL_ALPHABET *abc, ESL_MSA *rmsa, MSA_STAT *mrstat, ESL_MSA *emsa, MSA_STAT *mestat, float sc, 
		 float treeavgt, float time, int lcu_r, int lcu_e, char *errbuf, int verbose)
{
  double   sen    = 0.0;
  double   ppv    = 0.0;
  double   F      = 0.0;
  double   TC     = 0.0;
  double   SPE    = 0.0;
  int      alen;
  int      trueh  = 0;
  int      foundh = 0;
  int      tph    = 0;
  int      cac    = 0; // correctly aligned columns
  int      ac     = 0; // aligned columns in reference alignment
  int      nhr    = 0; // number of nonhomologous  positions in reference
  int      nhe    = 0; // number of nonhomologous  positions in inferred
  int      hr     = 0; // number of homologous  positions in reference
  int      he     = 0; // number of homologous  positions in inferred
  int      hre    = 0; // number of homologous  positions in reference and inferred
  int      status; 
  
  status = FastSP_Run(rmsa, emsa, &tph, &trueh, &foundh, &cac, &ac, &sen, &ppv, &F, &TC, lcu_r, lcu_e, errbuf, verbose);
  if (status != eslOK) goto ERROR;

  alen = (emsa)? (int)emsa->alen : 0.0;

  /* nonhomologies */
  status = msamanip_NonHomologous(abc, rmsa, emsa, &nhr, &nhe, &hr, &he, &hre, errbuf);
  if (status != eslOK) goto ERROR;

  SPE = (nhr > 0)? (float)nhe/(float)nhr : 0.0;
 
  if (method) fprintf(benchfp, "%s         %10d %10d %10d %10d %10d  %2.2f %2.2f %2.2f %2.2f    %f    %f    %s    %" PRId64 "  %d  %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d\n", 
		      msaname, tph, trueh, foundh, nhe, nhr, mrstat->avgid, mestat->avgid, mrstat->avgmatch, mestat->avgmatch, treeavgt, time, method, 
		      rmsa->alen, alen, mrstat->avgsqlen, mrstat->stdsqlen, mestat->avgsqlen, mestat->stdsqlen, 
		      mrstat->avginum, mrstat->stdinum, mestat->avginum, mestat->stdinum, mrstat->avgilen, mrstat->stdilen, mestat->avgilen, mestat->stdilen, sc, hr, he, hre);
  else        fprintf(benchfp, "%s         %10d %10d %10d %10d %10d   %2.2f %2.2f %2.2f %2.2f    %f    %f    NA    %" PRId64 "  %d  %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d\n",      
		      msaname, tph, trueh, foundh, nhe, nhr, mrstat->avgid, mestat->avgid, mrstat->avgmatch, mestat->avgmatch, treeavgt, time, 
		      rmsa->alen, alen, mrstat->avgsqlen, mrstat->stdsqlen, mestat->avgsqlen, mestat->stdsqlen,  
		      mrstat->avginum, mrstat->stdinum, mestat->avginum, mestat->stdinum, mrstat->avgilen, mrstat->stdilen, mestat->avgilen, mestat->stdilen, sc, hr, he, hre);

  if (1||verbose) fprintf(stdout, "#%s         %10d %10d %10d %10d %10d |  %.4f   %.4f   %.4f   %.4f | %2.2f %2.2f %2.2f %2.2f | %.4f     %.2f | %" PRId64 "  %d  %.2f +/- %.2f %.2f +/- %.2f %.4f | %d %d %d\n",      
			  msaname, tph, trueh, foundh, nhe, nhr, sen, ppv, F, SPE, mrstat->avgid, mestat->avgid, mrstat->avgmatch, mestat->avgmatch, 
			  treeavgt, time, 
			  rmsa->alen, alen, mrstat->avgsqlen, mrstat->stdsqlen, mestat->avgsqlen, mestat->stdsqlen, sc, hr, he, hre);
   return eslOK;

 ERROR:
  return status;
}


int
fastsp_parse(char *fastspfile, int *ret_tph, int *ret_trueh, int *ret_foundh, int *ret_cac, int *ret_ac, double *ret_sen, double *ret_ppv, double *ret_F, 
	     double *ret_TC, int flip, char *errbuf, int verbose)
{
  ESL_FILEPARSER *efp = NULL;
  char           *tok1;
  char           *tok;
  int             trueh  = 0;
  int             foundh = 0;
  int             tph    = 0;
  int             cac    = 0; // correctly aligned columns
  int             ac     = 0; // aligned columns in reference alignment
  double          sen;
  double          ppv;
  double          F;
  double          TC;
  int             aux;
  int             lnum = 0;
  int             status;

  if (esl_fileparser_Open(fastspfile, NULL, &efp) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "fastspfile open failed");
  esl_fileparser_SetCommentChar(efp, '#');
  
  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      /* get first token of each line */
      esl_fileparser_GetTokenOnLine(efp, &tok1, NULL);
      
      if (strcmp(tok1, "Number") == 0) {
	lnum ++;
	// need to get the last token of the line
	while((status = esl_fileparser_GetTokenOnLine(efp, &tok, NULL)) == eslOK) {
	  if (lnum == 1) tph    = atoi(tok); // Number of shared homologies: 309
	  if (lnum == 2) trueh  = atoi(tok); // Number of homologies in the reference alignment: 450
	  if (lnum == 3) foundh = atoi(tok); // Number of homologies in the estimated alignment: 464
	  if (lnum == 4) cac    = atoi(tok); // Number of correctly aligned columns: 51
	  if (lnum == 5) ac     = atoi(tok); // Number of aligned columns in ref. alignment: 193
	}
      }    
    }
  
  if (flip) {
    aux    = trueh;
    trueh  = foundh;
    foundh = aux;
  }
  sen = (trueh  > 0.)? (double)tph/(double)trueh  : 0.0;
  ppv = (foundh > 0.)? (double)tph/(double)foundh : 0.0;
  F = calculateF(sen, ppv);
  TC = (ac > 0)? (double)cac/(double)ac : 0.0;

  *ret_tph    = tph;
  *ret_trueh  = trueh;
  *ret_foundh = foundh;
  *ret_cac    = cac;
  *ret_ac     = ac;
  *ret_sen    = sen;
  *ret_ppv    = ppv;
  *ret_F      = F;
  *ret_TC     = TC;
  
  esl_fileparser_Close(efp);
  return eslOK;

 ERROR:
  if (efp) esl_fileparser_Close(efp);
  return status;
}


static double
calculateF(double sen, double ppv)
{
  double F;

  if (isnan(sen) || isnan(ppv)) F = 0.0;
  else if (sen+ppv == 0.0)      F = 0.0;
  else                          F = 2.0 * sen * ppv / (sen + ppv);

  return F;
}
