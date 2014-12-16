/*  orthologs - given a msa find ortholog sequences
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
#include "esl_cluster.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "msatree.h"
#include "orthologs.h"
#include "ssifile.h"

struct ortho_param_s {
  char     *ssiriotbl;
  double    minorthsc;
  char    **riotaxalist;
  int       riontaxa;
  int       verbose;
};

static int ortho_linkage(const void *v1, const void *v2, const void *p, int *ret_link);

ORTHO *
ortho_Create(int nc, char **taxalist, int ntaxa)
{
  ORTHO *ortho = NULL;
  int    n;
  int    c;
  int    status;

  /* Contract verification  */
  ESL_DASSERT1((nc >= 1));

  /* 1st allocation round  */
  ESL_ALLOC(ortho, sizeof(ORTHO) * nc);
  
  for (c = 0; c < nc; c ++) {
    ortho[c].useme = NULL;
    
    /* all orthologs untill shown otherwise */
    ESL_ALLOC(ortho[c].taxalist, sizeof(char *) * ntaxa);
    for (n = 0; n < ntaxa; n ++) ortho[c].taxalist[n] = NULL;
    if (taxalist) {
      for (n = 0; n < ntaxa; n ++) 
	esl_strdup(taxalist[n], -1, &(ortho[c].taxalist[n]));
    }
    
    /* all orthologs untill shown otherwise */
    ESL_ALLOC(ortho[c].useme, sizeof(int) * ntaxa);
    for (n = 0; n < ntaxa; n ++) ortho[c].useme[n] = TRUE;
    
    ortho[c].ntaxa     = ntaxa;
    ortho[c].no        = 0;
    ortho[c].orthodist = eslINFINITY;
  }
  return ortho;
  
 ERROR:
  ortho_Destroy(nc, ortho);
  return NULL;
}

void
ortho_Destroy(int nc, ORTHO *ortho)
{
  int c;
  int n;

  if (ortho == NULL) return;

  for (c = 0; c < nc; c ++) {
    if (ortho[c].useme) free(ortho[c].useme);
    if (ortho[c].taxalist) {
      for (n = 0; n < ortho[c].ntaxa; n ++) if (ortho[c].taxalist[n]) free(ortho[c].taxalist[n]);
      free(ortho[c].taxalist); 
    }
  }
  
  free(ortho);
  return;
}

int
ortho_FindClusters(char **sptreefile, FILE *outfp, ESL_RANDOMNESS *r, ESL_MSA *msa, char *outheader, ORTHO **ret_cortho, int *ret_nc, int nboot, int reboot, 
		   double minorthsc, char *errbuf, int verbose)
{
  char *treebootfile = NULL;
  int   status;

  /* (1) create the boostrapped gene trees */
  esl_sprintf(&treebootfile, "%s.treeboot", outheader);
  if ( !esl_FileExists(treebootfile) || reboot)
    if ((status = ortho_BoostrapTrees(r, msa, outheader, nboot, treebootfile, errbuf, verbose)) != eslOK) { printf("%s\n", errbuf); goto ERROR; }
 
  /* (2) run orthologs program */
  if ((status = ortho_RunRIO(sptreefile, outfp, msa, treebootfile, ret_cortho, ret_nc, minorthsc, errbuf, verbose)) != eslOK) 
    ESL_XFAIL(status, errbuf, "failed to run RIO");
  
  free(treebootfile);
  return eslOK;
  
 ERROR:
  if (treebootfile) free(treebootfile);
  return status;
}

int
ortho_BoostrapTrees(ESL_RANDOMNESS *r, ESL_MSA *msa, char *outheader, int nboot, char *treebootfile, char *errbuf, int verbose)
{
  char           *msabootfile = NULL;
  FILE           *treebootfp = NULL;
  ESLX_MSAFILE   *msafp = NULL;
  ESL_TREE      **Tlist = NULL;
  ESL_MSA        *bmsa = NULL;
  int             n = 0;
  int             b;
  int             s;
  int             hstatus = eslOK;   
  int             status;

  /* create the boostrapped msas using SEQBOOT 
   * returns bootstrapped msas in stockholm format.
   */
  msabootfile = ortho_RunSEQBOOT(r, msa, outheader, nboot, errbuf, verbose);
  if (msabootfile == NULL) { status = eslEMEM; goto ERROR; }
  
  /* allocate for the boostrapped trees Tlist */
  ESL_ALLOC(Tlist, sizeof(ESL_TREE *) * nboot);
  for (b = 0; b < nboot; b++) Tlist[b] = NULL;
  
  /* Open the msaboot file */
  status = eslx_msafile_Open(&(msa->abc), msabootfile, NULL, eslMSAFILE_STOCKHOLM, NULL, &msafp);
  if (status != eslOK) eslx_msafile_OpenFailure(msafp, status);
  while ((hstatus = eslx_msafile_Read(msafp, &bmsa)) == eslOK) 
    {
      /* bit of hackery to restore seqnames.
       * seqboot shortens them to 10 characters. 
       * The order of seqs has not change by running seqboot.
       * We want to restore the original names so those can be passed to the tree.
       */
      for (s = 0; s < msa->nseq; s++) {
	free(bmsa->sqname[s]);
	esl_strdup(msa->sqname[s], -1, &(bmsa->sqname[s]));
      }
      
      /* create corresponding T.
       * Don't need to do anything about the rooting of the tree. 
       * RIO does re-roots each gene tree as to minimize the sum of duplications */
      if ((status = Tree_CalculateExtFromMSA(bmsa, &(Tlist[n]), FALSE, errbuf, FALSE)) != eslOK) { printf("%s\n", errbuf); goto ERROR; }
      n ++;
      esl_msa_Destroy(bmsa); bmsa = NULL;
  }  
  if (n != nboot) { printf("wrong number of boostrapped msas %d, should be %d\n", n, nboot); status = eslFAIL; goto ERROR; }
  eslx_msafile_Close(msafp); 

  /* write the treeboot file */
  remove(treebootfile);
  for (b = 0; b < nboot; b++) {
    if ((treebootfp = fopen(treebootfile, "a+")) == NULL) { printf("failed to open %s\n", treebootfile); goto ERROR; }
    esl_tree_WriteNewick(treebootfp, Tlist[b]);
    fclose(treebootfp);
  }

  remove(msabootfile); free(msabootfile);
 if (Tlist) {
    for (b = 0; b < nboot; b++) if (Tlist[b]) esl_tree_Destroy(Tlist[b]);
    free(Tlist);  
  }
  return eslOK;

 ERROR:
  if (Tlist) {
    for (b = 0; b < nboot; b++) if (Tlist[b]) esl_tree_Destroy(Tlist[b]);
    free(Tlist);
  }
  if (bmsa) esl_msa_Destroy(bmsa); 
  remove(msabootfile);
   if (msabootfile)  free(msabootfile);
  if (treebootfile) free(treebootfile);
  if (Tlist) {
    for (b = 0; b < nboot; b++) if (Tlist[b]) esl_tree_Destroy(Tlist[b]);
    free(Tlist);  
  }

  return status;
}

char *
ortho_RunSEQBOOT(ESL_RANDOMNESS *r, ESL_MSA *msa, char *outheader, int N, char *errbuf, int verbose)
{
  char  *msabootfile = NULL;  
  char  *phylipcmd = NULL;                 /* phylip command line file */
  char  *phylipin  = NULL;                 /* phylip infile */
  char  *phylipout = NULL;                 /* phylip outfile */
  char  *s;
  char  *args = NULL;
  FILE  *msafp = NULL;
  FILE  *cmdfp = NULL;
  int    oddran;
  int    status;
  
  if ("EVOHMMDIR" == NULL)               { status = eslENOTFOUND; goto ERROR; }
  if ((s = getenv("EVOHMMDIR")) == NULL) { status = eslENOTFOUND; goto ERROR; }

 /* msaboot  file */
  esl_sprintf(&msabootfile, "%s.boot",      outheader);
  esl_sprintf(&phylipcmd,   "%s.phylipcmd", msabootfile);
  esl_sprintf(&phylipin,    "%s.phylipin",  msabootfile);
  esl_sprintf(&phylipout,   "outfile");
  
  /* save msa in phylip format */
  if ((msafp = fopen(phylipin, "w")) == NULL) { printf("failed to open %s for writting\n", phylipin); goto ERROR; }
  if ((status = eslx_msafile_Write(msafp, msa, eslMSAFILE_PHYLIP)) != eslOK) goto ERROR; 
  fclose(msafp);

  oddran = (int)(esl_random(r) * 10000.);
  if (oddran%2 == 0) oddran ++;
  /* create command line file */
  if ((cmdfp = fopen(phylipcmd, "w")) == NULL) { printf("failed to open %s for writting\n", phylipcmd); goto ERROR; }
  fprintf(cmdfp, "%s\n", phylipin);  /* the input file name */
  fprintf(cmdfp, "2\n");             /* don't print details of run */
  fprintf(cmdfp, "R\n");             /* number of samples   */
  fprintf(cmdfp, "%d\n", N);         /* number of samples   */
  fprintf(cmdfp, "y\n");
  fprintf(cmdfp, "%d\n", oddran);    /* odd number to seed a random number generator */
  fprintf(cmdfp, "y\n");
  fclose(cmdfp);
  
  /* run seqboot */
  esl_sprintf(&args, "%s/lib/phylip-3.69/src/seqboot < %s > foo\n", s, phylipcmd);
  system(args);
  free(args); args = NULL;
  
  /* convert boostrapped msa (outfile) to stockholm format */
  esl_sprintf(&args, "%s/lib/hmmer4/lib/easel/miniapps/esl-reformat stockholm %s > %s\n", s, phylipout, msabootfile);
  system(args);

  remove(phylipcmd);
  remove(phylipin);
  remove(phylipout);
  
  free(args);
  free(phylipcmd);
  free(phylipin);
  free(phylipout);
  
  return msabootfile;
  
 ERROR:
  remove(phylipin);
  remove(phylipout);
  remove(phylipcmd);
  if (args)      free(args);
  if (phylipcmd) free(phylipcmd);
  if (phylipin)  free(phylipin);
  if (phylipout) free(phylipout);
   return NULL;
}


int
ortho_RunRIO(char **sptreefile, FILE *outfp, ESL_MSA *msa, char *treebootfile, ORTHO **ret_cortho, int *ret_nc, double minorthsc, char *errbuf, int verbose)
{
  char      *riotbl = NULL;
  char      *rioout = NULL;
  char      *s = NULL;
  char      *args = NULL;
  int        forester_version = 1019;
  int        rootatmidpoint = TRUE;
  int        status;
  
  if ("EVOHMMDIR" == NULL)               { status = eslENOTFOUND; goto ERROR; }
  if ((s = getenv("EVOHMMDIR")) == NULL) { status = eslENOTFOUND; goto ERROR; }

  if (sptreefile == NULL)                { status = eslENOTFOUND; goto ERROR; }
 
#if 0
  /* TODO: our own species tree -- use Christian's for now */
  if (*sptreefile) free(*sptreefile);
  esl_sprintf(sptreefile, "%s/lib/RIO/species.xml", s);
#endif

  esl_sprintf(&riotbl, "%s.tbl", treebootfile);
  esl_sprintf(&rioout, "%s.out", treebootfile);
  
  /* run RIO */
  if (rootatmidpoint) {
    if (verbose) fprintf(outfp, "\njava -Xmx2048m -cp %s/lib/RIO/forester_%d.jar org.forester.application.rio -r=midpoint %s %s %s %s\n", s, forester_version, treebootfile, *sptreefile, riotbl, rioout);
    esl_sprintf(&args,            "java -Xmx2048m -cp %s/lib/RIO/forester_%d.jar org.forester.application.rio -r=midpoint %s %s %s %s>foo\n", s, forester_version, treebootfile, *sptreefile, riotbl, rioout);
  }
  else {
    if (verbose) fprintf(outfp, "\njava -Xmx2048m -cp %s/lib/RIO/forester_%d.jar org.forester.application.rio %s %s %s %s\n", s, forester_version, treebootfile, *sptreefile, riotbl, rioout);
    esl_sprintf(&args,            "java -Xmx2048m -cp %s/lib/RIO/forester_%d.jar org.forester.application.rio %s %s %s %s>foo\n", s, forester_version, treebootfile, *sptreefile, riotbl, rioout);
  }
  system(args);
  free(args); args = NULL;

  /* get orthologs from riotbl */
  if ((status = ortho_LinkageRIOTable(riotbl, outfp,  msa, ret_cortho, ret_nc, minorthsc, errbuf, verbose) ) != eslOK) ESL_XFAIL(status, errbuf, "failed to find orthologs");

  remove(riotbl);
  remove(rioout);
  free(riotbl);
  free(rioout);
  return eslOK;

 ERROR:
  if (args) free(args);
  remove(riotbl);
  remove(rioout);
  if (riotbl) free(riotbl);
  if (rioout) free(rioout);
  return status;
}

int
ortho_LinkageRIOTable(char *riotbl, FILE *outfp, ESL_MSA *msa, ORTHO **ret_cortho, int *ret_nc, double minorthsc, char *errbuf, int verbose)
{
  char   *ssiriotbl = NULL;
  char  **riotaxalist = NULL;
  FILE   *fp = NULL;
  int     riontaxa;                     /* sequences actually analyzed by RIO */
  int     x;
  int     status;

  if (riotbl == NULL) { return eslOK; }

 /* Create an ssi index of the riotbl */
  status = ssi_riotbl_Create(riotbl, &ssiriotbl, &riotaxalist, &riontaxa);
  if (status != eslOK) { printf("ssi_riotbl_Create failed for file %s\n", riotbl); status = eslFAIL; goto ERROR; }
  if (verbose) fprintf(outfp, "ssifile %s created\n", ssiriotbl);
  
  fprintf(outfp, "\nRIO: analizes %d / %d sequences\n", riontaxa, msa->nseq);
  fprintf(outfp, "RIO: cutoff for orthology %f\n", minorthsc);

  if ((fp = fopen(riotbl, "r")) == NULL) { printf("failed to open %s\n", riotbl); status = eslFAIL; goto ERROR; }

#if 1
  if (ortho_CompleteLinkage(ssiriotbl, outfp, msa, riotaxalist, riontaxa, ret_cortho, ret_nc, minorthsc, errbuf, verbose) != eslOK)
  { printf("failed to do complete-linkage clustering\n"); status = eslFAIL; goto ERROR; }
#else
  if (ortho_SingleLinkage(ssiriotbl, outfp, msa, riotaxalist, riontaxa, ret_cortho, ret_nc, minorthsc, errbuf, verbose) != eslOK)
  { printf("failed to do single-linkage clustering\n"); status = eslFAIL; goto ERROR; }
#endif

  fclose(fp);
  remove(ssiriotbl); free(ssiriotbl);
  for (x = 0; x < riontaxa; x ++) free(riotaxalist[x]);
  free(riotaxalist);
  return eslOK;

 ERROR:
  remove(ssiriotbl);
  if (riotaxalist) {
    for (x = 0; x < riontaxa; x ++) if (riotaxalist[x]) free(riotaxalist[x]);
    free(riotaxalist);
  }
  return status;
}
 
int
ortho_CompleteLinkage(char *ssiriotbl, FILE *outfp, ESL_MSA *msa, char **riotaxalist, int riontaxa, ORTHO **ret_cortho, int *ret_nc, double minorthsc, char *errbuf, int verbose)
{
  ESL_DMATRIX *O = NULL;
  ESL_TREE    *LT = NULL; /* linkage tree */
  float        orthsc;
  int          n, m;
  int          status;

  O = esl_dmatrix_Create(msa->nseq, msa->nseq);
  if (O == NULL) { status = eslFAIL; goto ERROR; }

  for (n = 0; n < msa->nseq-1; n ++) {
    O->mx[n][n] = 0.0;
    for (m = n+1; m < msa->nseq; m ++) {
      ssi_riotbl_GetOrthsc(ssiriotbl, riotaxalist, riontaxa, msa->sqname[n], msa->sqname[m], &orthsc, verbose);
      O->mx[n][m] = -log(orthsc);
      O->mx[m][n] = -log(orthsc);
    }
  }
  O->mx[msa->nseq-1][msa->nseq-1] = 0.0;

  if (verbose) esl_dmatrix_Dump(stdout, O, NULL, NULL);

  if (esl_tree_CompleteLinkage(O, &LT) != eslOK) { status = eslFAIL; goto ERROR; }
  if (verbose) Tree_Dump(stdout, LT, "CompleteLinkage Tree");

  if (ortho_LinkageTree2Clusters(outfp, msa, LT, ret_cortho, ret_nc, minorthsc, errbuf, verbose) != eslOK) { status = eslFAIL; goto ERROR; }

  /* cleanup and return */
  esl_dmatrix_Destroy(O);
  esl_tree_Destroy(LT);
  return eslOK;

 ERROR:
  if (O) esl_dmatrix_Destroy(O);
  if (LT) esl_tree_Destroy(LT);
  return status;
}

int
ortho_SingleLinkage(char *ssiriotbl, FILE *outfp, ESL_MSA *msa, char **riotaxalist, int riontaxa, ORTHO **ret_cortho, int *ret_nc, double minorthsc, char *errbuf, int verbose)
{
  ORTHO  *cortho = NULL;
  int    *workspace  = NULL;
  int    *assignment = NULL;
  int     nc;
  int     ncs = 0;
  int     c;
  int     i;
  int     u;
  int     status;
  struct ortho_param_s param;

  if (ssiriotbl == NULL) return eslOK;

  /* Allocations */
  ESL_ALLOC(workspace,  sizeof(int) * msa->nseq * 2);
  ESL_ALLOC(assignment, sizeof(int) * msa->nseq);

  /* call to SLC API: */
  param.minorthsc   = minorthsc;
  param.riotaxalist = riotaxalist;
  param.riontaxa    = riontaxa;
  param.ssiriotbl   = ssiriotbl;
  param.verbose     = verbose;
  status = esl_cluster_SingleLinkage((void *) msa->sqname, (size_t) msa->nseq, sizeof(char *),
				     ortho_linkage, (void *) &param, 
				     workspace, assignment, &nc);

  cortho = ortho_Create(nc, msa->sqname, msa->nseq);
  if (cortho == NULL) { status = eslEMEM; goto ERROR; }

  for (i = 0; i < msa->nseq; i++)
    cortho[assignment[i]].no ++;
     
   for (c = 0; c < nc; c ++) {
     if (cortho[c].no > 1) {
	fprintf(outfp, "cluster[%d] = %d\n", c, cortho[c].no);
	ncs ++;
     }
   }
   fprintf(outfp, "ORTHOLOGY: Found %d single-linkage clusters with more than 1 seq\n", ncs);
 
  for (c = 0; c < nc; c ++) { 
    u = 0;
    if (cortho[c].no > cortho[c].ntaxa) { printf("found gibberish numbers of ortologs\n"); status = eslFAIL; goto ERROR; }
    if (cortho[c].no == 0)              { printf("failed to find any orthologs\n");        status = eslFAIL; goto ERROR; }
    
    for (i = 0; i < cortho[c].ntaxa; i++) cortho[c].useme[i] = (assignment[i] == c) ? TRUE : FALSE;
    for (i = 0; i < cortho[c].ntaxa; i++) { if (cortho[c].useme[i]) u ++; }
    if (u != cortho[c].no) { printf("failed creating ortho list has %d, should be %d\n", u, cortho[c].no); status = eslFAIL; goto ERROR; }
  }

  *ret_cortho = cortho;
  *ret_nc     = nc;
  
  /* cleanup and return */
  free(workspace);
  free(assignment);
  return eslOK;

 ERROR:
  if (workspace  != NULL) free(workspace);
  if (assignment != NULL) free(assignment);
   return status;
}

int
ortho_LinkageTree2Clusters(FILE *outfp, ESL_MSA *msa, ESL_TREE *T, ORTHO **ret_cortho, int *ret_nc, double minorthsc, char *errbuf, int verbose)
{
  ESL_STACK   *vs = NULL;
  ESL_STACK   *cs = NULL;
  ORTHO       *cortho = NULL;
  double      *odist = NULL;
  int         *assignment = NULL;
  double       minsc = -log(minorthsc);
  int          nc = 0;
  int          ncs = 0;
  int          no;
  int          v;
  int          vv;
  int          i;
  int          c;
  int          u;
  int          status;

  if ((status = esl_tree_SetCladesizes(T))  != eslOK) ESL_XFAIL(status, errbuf, "failed to SetCladesize");
  
  /* Allocations */
  ESL_ALLOC(assignment, sizeof(int) * T->N);
  esl_vec_ISet(assignment, T->N, -1);

  /* create a stack, and put the root in the stack */
  if ((vs = esl_stack_ICreate()) == NULL) { status = eslEMEM; goto ERROR; };
  if (esl_stack_IPush(vs, 0) != eslOK)     { status = eslEMEM; goto ERROR; };

  while (esl_stack_IPop(vs, &v) == eslOK) 
    {                                    
      /* this is a linkage tree so left/right branch values are the same */
      if (T->ld[v] <= minsc) { /* a cluster of orthologs */

	/* make the assignment for all sequences under this node */
	if ((cs = esl_stack_ICreate()) == NULL) { status = eslEMEM; goto ERROR; };
	if (esl_stack_IPush(cs, v) != eslOK)    { status = eslEMEM; goto ERROR; };
	while (esl_stack_IPop(cs, &vv) == eslOK) 
	  {
	    if (T->left[vv]  > 0) esl_stack_IPush(cs, T->left[vv]);  else assignment[-T->left[vv]]  = nc;
	    if (T->right[vv] > 0) esl_stack_IPush(cs, T->right[vv]); else assignment[-T->right[vv]] = nc;
	  }
	esl_stack_Destroy(cs); cs = NULL;
	
	if (odist == NULL) ESL_ALLOC  (odist,  sizeof(double) * (nc+1));
	else               ESL_REALLOC(odist,  sizeof(double) * (nc+1));
	odist[nc] = T->ld[v];
	nc ++;
      }
      else {
 	if (T->left[v]  > 0 && esl_stack_IPush(vs, T->left[v])  != eslOK) { status = eslEMEM; goto ERROR; }; 
 	if (T->right[v] > 0 && esl_stack_IPush(vs, T->right[v]) != eslOK) { status = eslEMEM; goto ERROR; }; 
      }
   }

  cortho = ortho_Create(nc, msa->sqname, msa->nseq);
  if (cortho == NULL) { status = eslEMEM; goto ERROR; }

  for (i = 0; i < msa->nseq; i++)
    if (assignment[i] >= 0) cortho[assignment[i]].no ++;
  
   for (c = 0; c < nc; c ++) 
     if (cortho[c].no > 1) 
	ncs ++;
   fprintf(outfp, "ORTHOLOGY: Found %d clusters with more than 1 seq\n", ncs);

  for (c = 0; c < nc; c ++) { 
    u = 0;
    if (cortho[c].no > cortho[c].ntaxa) { printf("found gibberish numbers of ortologs\n"); status = eslFAIL; goto ERROR; }
    if (cortho[c].no == 0)              { printf("failed to find any orthologs\n");        status = eslFAIL; goto ERROR; }

    cortho[c].orthodist = odist[c];
    
    for (i = 0; i < cortho[c].ntaxa; i++) cortho[c].useme[i] = (assignment[i] == c) ? TRUE : FALSE;
    for (i = 0; i < cortho[c].ntaxa; i++) { if (cortho[c].useme[i]) u ++; }
    if (u != cortho[c].no) { printf("failed creating ortho list has %d, should be %d\n", u, cortho[c].no); status = eslFAIL; goto ERROR; }
  }
  
  for (c = 0; c < nc; c ++) {
    if (cortho[c].no > 1) {
      ncs ++;

      if (1||verbose) 
	{
	  no = 0;
	  fprintf(outfp, "cluster[%d] = %d || orthosc %f\n", c, cortho[c].no, exp(-cortho[c].orthodist));
#if 0
	  for (i = 0; i < cortho[c].ntaxa; i ++) {
	    if (cortho[c].useme[i]) {
	      fprintf(outfp, "oseq[%d]: %s\n", no++, cortho[c].taxalist[i]);
	    }
	  }
#endif
	}
    }
  }
  
  
  *ret_cortho = cortho;
  *ret_nc     = nc;

  free(assignment);
  free(odist);
  esl_stack_Destroy(vs);
  return eslOK;

 ERROR:
  if (assignment) free(assignment);
  if (odist) free(odist);
  ortho_Destroy(nc, cortho);
  if (vs) esl_stack_Destroy(vs);
  if (cs) esl_stack_Destroy(cs);
  return status;
}


/*****************************************************************
 * 2. Internal functions, interface to the clustering API
 *****************************************************************/

/* Definition of linkage  (rioboostrap >= minorthsc): */
static int
ortho_linkage(const void *v1, const void *v2, const void *p, int *ret_link)
{
  char  *tx1 = *(char **) v1;
  char  *tx2 = *(char **) v2;
  struct ortho_param_s *param = (struct ortho_param_s *) p;
  float  orthsc;
  int    status = eslOK;

  ssi_riotbl_GetOrthsc(param->ssiriotbl, param->riotaxalist, param->riontaxa, tx1, tx2, &orthsc, param->verbose);

  *ret_link = (orthsc >= param->minorthsc ? TRUE : FALSE); 
  return status;
}
  
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
