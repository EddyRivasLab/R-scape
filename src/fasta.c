/*  fasta -- function to integrate FASTA msa output
 *
 * ER, Fri Oct 31 15:55:28 EDT 2014 [Janelia] 
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
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "msamanip.h"
#include "minimize.h"
#include "miscellaneous.h"
#include "fasta.h"


int
SSEARCH_Align(char *version, const ESL_MSA *msa, ESL_MSA **ret_ssearchmsa, float *ret_sc, char *mx, int fasta_f, int fasta_g, char *errbuf, int verbose)
{
  ESL_MSAFILE    *afp = NULL;
  FILE            *fp = NULL;
  ESL_MSA         *ssearchmsa = NULL;
  ESL_SQ          *sq = NULL;
  char             tmpdbfile[16]    = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpqueryfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char             tmpmsafile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char            *args = NULL;
  char            *ssearchopts = NULL;
  char            *s = NULL;
  float            sc = -eslINFINITY;
  int              bestd = 0;
  int              d;
  int              status;
  
  if (msa->nseq != 2)  ESL_XFAIL(status, errbuf, "SSEARCH_align() only works for pair of sequences for now");
  
  /* db in FASTA format */
  if ((status = esl_tmpfile_named(tmpdbfile,  &fp))   != eslOK) ESL_XFAIL(status, errbuf, "failed to create db file");
  esl_sq_FetchFromMSA(msa, 0, &sq);
  if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write db  file\n");
  fclose(fp); 
  esl_sq_Destroy(sq); sq = NULL;
  
  /* query in FASTA format */
  if ((status = esl_tmpfile_named(tmpqueryfile,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create query file");
  esl_sq_FetchFromMSA(msa, 1, &sq);
  if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write query file\n");
  fclose(fp);
  esl_sq_Destroy(sq); sq = NULL;
  
  /* run SSEARCH */
  if ("FASTADIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("FASTADIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpoutfile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  fclose(fp);

  //ssearch36 options
  // -q       quiet
  // -E 10000 evalue
  // -b 1     to report only the highest score
  // -d 1     number of alignments shown (limited by -E by default)
  // 
  if (strcmp(mx, "default") == 0) esl_sprintf(&ssearchopts, "-q -b 1 -d 1 ");
  else {
    if (fasta_f == 0 || fasta_g == 0) ESL_XFAIL(status, errbuf, "need to provide customized gap-open and gap-extent penalties\n");
    esl_sprintf(&ssearchopts, "-q -b 1 -d 1 -s %s -f %d -g %d ", mx, fasta_f, fasta_g);
  }

  esl_sprintf(&args, "%s/%s/bin/ssearch36 %s %s %s > %s", s, version, ssearchopts, tmpqueryfile, tmpdbfile, tmpoutfile);
  system(args);
  printf("%s\n", args);

  /* parse ssearch output to get alignment */ 
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpmsafile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create msa file");
  fclose(fp);
  // report best alignment of best hit (-a) and add unaligned flanking regions (-f)
  esl_sprintf(&args, "%s/benchmarks/scripts/fasta2afa.pl -a -f -q %s -t %s  %s %s ", s, tmpqueryfile, tmpdbfile, tmpoutfile, tmpmsafile);
  system(args);
  printf("%s\n", args);

  /* convert to msa */
  if (esl_msafile_Open(NULL, tmpmsafile, NULL, eslMSAFILE_AFA, NULL, &afp) != eslOK) ESL_XFAIL(status, errbuf, "Failed to open AFA file\n");
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);
  if (esl_msafile_Read(afp, &ssearchmsa) != eslOK) ESL_XFAIL(status, errbuf, "Failed to read AFA file\n");
  if (status != eslOK) esl_msafile_ReadFailure(afp, status);
  esl_msafile_Close(afp);

  *ret_sc = sc;
  if (ret_ssearchmsa) *ret_ssearchmsa = ssearchmsa;
 
  //remove(tmpdbfile);
  //remove(tmpqueryfile);
  //remove(tmpoutfile);
  //remove(tmpmsafile);
  
  if (sq) esl_sq_Destroy(sq);
  if (ssearchopts) free(ssearchopts);
  if (args) free(args);
  return eslOK;
  
 ERROR:
  //remove(tmpdbfile);
  //remove(tmpqueryfile);
  //remove(tmpoutfile);
  //remove(tmpmsafile);
  if (sq) esl_sq_Destroy(sq);
  if (ssearchopts) free(ssearchopts);
  if (args) free(args);
  return status;  
}


 
