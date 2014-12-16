/*  blast -- function to integrate BLAST msa output
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
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "ncbiblast.h"

int
NCBIBLAST_Align(const ESL_MSA *msa, int wordsize, ESL_MSA **ret_blastmsa, char *errbuf, int verbose)
{
  ESLX_MSAFILE *afp = NULL;
  FILE         *fp = NULL;
  ESL_MSA      *blastmsa = NULL;
  ESL_SQ       *sq = NULL;
  char          tmpdbfile[16]    = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpqueryfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpoutfile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char          tmpmsafile[16]   = "esltmpXXXXXX"; /* tmpfile template */
  char         *tmpdbfilephr = NULL; /* tmpfile template */
  char         *tmpdbfilepin = NULL; /* tmpfile template */
  char         *tmpdbfilepsq = NULL; /* tmpfile template */
  char         *args = NULL;
  char         *blastopts = NULL;
  char         *s = NULL;
  int           status;
  
  if (msa->nseq != 2)  ESL_XFAIL(status, errbuf, "NCBIBLAS_align() only works for pair of sequences for now");
  
  /* BLAST db in FASTA format */
  if ((status = esl_tmpfile_named(tmpdbfile,  &fp))   != eslOK) ESL_XFAIL(status, errbuf, "failed to create db file");
  esl_sq_FetchFromMSA(msa, 0, &sq);
  if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write db  file\n");
  fclose(fp); 
  esl_sq_Destroy(sq); sq = NULL;
  
  /* BLAST query in FASTA format */
  if ((status = esl_tmpfile_named(tmpqueryfile,  &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create query file");
  esl_sq_FetchFromMSA(msa, 1, &sq);
  if ((status = esl_sqio_Write(fp, sq, eslSQFILE_FASTA, FALSE)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write query file\n");
  fclose(fp);
  esl_sq_Destroy(sq); sq = NULL;
  
  /* index the db */
  if ("NCBIBLASTDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("NCBIBLASTDIR")) == NULL) return eslENOTFOUND;
  esl_sprintf(&args, "%s/makeblastdb -dbtype prot -in %s > /dev/null", s, tmpdbfile);
  system(args);
  esl_sprintf(&tmpdbfilephr, "%s.phr", tmpdbfile);
  esl_sprintf(&tmpdbfilepin, "%s.pin", tmpdbfile);
  esl_sprintf(&tmpdbfilepsq, "%s.psq", tmpdbfile);

  /* run NCBIBLAST */
  if ((status = esl_tmpfile_named(tmpoutfile, &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  if (wordsize > 0)  esl_sprintf(&blastopts, "-num_threads 1 -num_descriptions 9999 -evalue 10000 -word_size %d ", wordsize);
  else               esl_sprintf(&blastopts, "-num_threads 1 -num_descriptions 9999 -evalue 10000 ");
  esl_sprintf(&args, "%s/blastp -db %s -query %s %s > %s", s, tmpdbfile, tmpqueryfile, blastopts, tmpoutfile);
  system(args);
  printf("%s\n",args);
  fclose(fp);
  
  /* parse ncbiblast output to get alignment */ 
  if ("EVOHMMDIR" == NULL)               return eslENOTFOUND;
  if ((s = getenv("EVOHMMDIR")) == NULL) return eslENOTFOUND;
  if ((status = esl_tmpfile_named(tmpmsafile, &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create msa file");
  fclose(fp);
  // report best alignment of best hit (-a) and add unaligned flanking regions (-f)
  esl_sprintf(&args, "%s/benchmarks/scripts/blast2afa.pl -a -f -q %s -t %s  %s %s ", s, tmpqueryfile, tmpdbfile, tmpoutfile, tmpmsafile);
  system(args);

  /* convert to msa */
  if (eslx_msafile_Open(NULL, tmpmsafile, NULL, eslMSAFILE_AFA, NULL, &afp) != eslOK) ESL_XFAIL(status, errbuf, "Failed to open AFA file\n");
  if (status != eslOK) eslx_msafile_OpenFailure(afp, status);
  if (eslx_msafile_Read(afp, &blastmsa) != eslOK) ESL_XFAIL(status, errbuf, "Failed to read AFA file\n");
  if (status != eslOK) eslx_msafile_ReadFailure(afp, status);
  eslx_msafile_Close(afp);
  
  remove(tmpdbfile);
  remove(tmpdbfilephr);
  remove(tmpdbfilepin);
  remove(tmpdbfilepsq);
  remove(tmpqueryfile);
  remove(tmpoutfile);
  remove(tmpmsafile);
  
  *ret_blastmsa = blastmsa;

  if (sq) esl_sq_Destroy(sq);
  if (blastopts) free(blastopts);
  if (args) free(args);
  free(tmpdbfilephr);
  free(tmpdbfilepin);
  free(tmpdbfilepsq);
  return eslOK;
  
 ERROR:
  remove(tmpdbfile);
  remove(tmpdbfilephr);
  remove(tmpdbfilepin);
  remove(tmpdbfilepsq);
  remove(tmpqueryfile);
  remove(tmpoutfile);
  remove(tmpmsafile);
  if (sq) esl_sq_Destroy(sq);
  if (blastopts) free(blastopts);
  if (args) free(args);
  if (tmpdbfilephr) free(tmpdbfilephr);
  if (tmpdbfilepin) free(tmpdbfilepin);
  if (tmpdbfilepsq) free(tmpdbfilepsq);
  return status;  
}
