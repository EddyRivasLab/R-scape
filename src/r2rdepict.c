/* r2rdepict.c */

#include "rscape_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "covariation.h"
#include "correlators.h"

static int r2r_depict_pdf(char *r2rfile, int verbose, char *errbuf);
static int r2r_depict_svg(char *r2rfile, int verbose, char *errbuf);
static int r2r_pseudoknot(ESL_MSA *msa, int nct, int **ctlist);

int
r2r_Depict(char *r2rfile, int r2rall, ESL_MSA *msa, int nct, int **ctlist, HITLIST *hitlist, int makepdf, int makesvg, int verbose, char *errbuf)
 {
  ESL_MSAFILE *afp = NULL;
  FILE        *fp = NULL;
  int         *ct;
  ESL_MSA     *r2rmsa = NULL;
  char         tmpinfile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char         tmpoutfile[16] = "esltmpXXXXXX"; /* tmpfile template */
  char         covtag[12]     = "cov_SS_cons";
  char        *covtag1;
  char        *args = NULL;
  char        *cmd  = NULL;
  char        *ssstr      = NULL;
  char        *covstr     = NULL;
  char        *prv_covstr = NULL;
  char        *tok = NULL;
  char        *tag = NULL;
  int          found;
  int          i;
  int          h;
  int          ih, jh;
  int          tagidx, tagidx2;
  int          idx;
  int          do_r2rcovmarkup = FALSE;
  int          s;
  int          status;

  if (r2rfile == NULL) return eslOK;

  ESL_ALLOC(ssstr, sizeof(char) * (msa->alen+1));
  for (s = 0; s < nct; s++) {
    ct = ctlist[s];
    /* first modify the ss to a simple <> format. R2R cannot deal with fullwuss 
     */
    esl_ct2simplewuss(ct, msa->alen, ssstr);

    // if the primary structure (s==0), replace the 'SS_cons' GC line with the new ss
    //
    // all the other substructures add as GC tag as
    //
    // #GC SS_cons_1 
    // #GC SS_cons_2
    //
    if (s == 0) strcpy(msa->ss_cons, ssstr);
    else {
      esl_sprintf(&tag, "SS_cons_%d", s);
      esl_msa_AppendGC(msa, tag, ssstr);
    }
  }
 
  /* R2R input and output in PFAM format (STOCKHOLM in one single block) */
  if ((status = esl_tmpfile_named(tmpinfile,  &fp))                     != eslOK) ESL_XFAIL(status, errbuf, "failed to create input file");
  if ((status = esl_msafile_Write(fp, (ESL_MSA *)msa, eslMSAFILE_PFAM)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to write PFAM file\n");
  fclose(fp);
  
  /* run R2R */
  if ((status = esl_tmpfile_named(tmpoutfile, &fp)) != eslOK) ESL_XFAIL(status, errbuf, "failed to create output file");
  fclose(fp);
  
  if (RSCAPE_BIN)         // look for the installed executable
    esl_sprintf(&cmd, "%s/r2r", RSCAPE_BIN);  
  else
    ESL_XFAIL(status, errbuf, "Failed to find R2R executable\n");
  
  esl_sprintf(&args, "%s --GSC-weighted-consensus %s %s 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1", cmd, tmpinfile, tmpoutfile);
  status = system(args);
  if (status == -1) ESL_XFAIL(status, errbuf, "Failed to run R2R\n");
  
  /* convert output to r2rmsa */
  if (esl_msafile_Open(NULL, tmpoutfile, NULL, eslMSAFILE_PFAM, NULL, &afp) != eslOK) esl_msafile_OpenFailure(afp, status);
  afp->format = eslMSAFILE_PFAM;
  if (esl_msafile_Read(afp, &r2rmsa) != eslOK) esl_msafile_ReadFailure(afp, status);
  esl_msafile_Close(afp);

  /* make a cov_cons_ss (prv_covstr) line according to our hitlist */
  if (msa->alen != r2rmsa->alen) ESL_XFAIL(eslFAIL, errbuf, "r2r has modified the alignment");
  for (i = 1; i <= msa->alen; i ++) esl_strcat(&prv_covstr, -1, ".", 1);
  
  if (hitlist) {
    for (i = 1; i <= msa->alen; i ++) {
      found = FALSE;
      for (h = 0; h < hitlist->nhit; h ++) {
	ih = hitlist->hit[h].i+1;
	jh = hitlist->hit[h].j+1;
	
	if ((i == ih || i == jh) && hitlist->hit[h].bptype == WWc) { 
	  esl_sprintf(&tok, "2"); 
	  found = TRUE; 
	}
	
	if (found) break;
      }
      if (!found) esl_sprintf(&tok, ".");
      
      if (i == 1) esl_sprintf(&covstr, "%s", tok);
      else        esl_sprintf(&covstr, "%s%s", prv_covstr, tok);
      
      free(prv_covstr); prv_covstr = NULL;
      esl_sprintf(&prv_covstr, "%s", covstr);
      free(tok); tok = NULL;
      free(covstr); covstr = NULL;
    }
  }
  
  /* add line #=GF R2R keep allpairs 
   * so that it does not truncate ss.
   * cannot use the standard esl_msa_AddGF:
   *             esl_msa_AddGF(msa, "R2R", -1, "keep allpairs", -1);
   * since it does not parse with r2r
   *
   * turns out the above solution can only deal with the  <> annotation
   */
  if (r2rall) {
    for (tagidx = 0; tagidx < r2rmsa->ngf; tagidx++) {
      esl_strchop(r2rmsa->gf[tagidx], -1);
      if (strcmp(r2rmsa->gf[tagidx], "keep all") == 0) break;
    }
    
    if (tagidx < r2rmsa->ngf) { //remove 
      for (idx = tagidx; idx < r2rmsa->ngf-1; idx++) {
	esl_sprintf(&r2rmsa->gf_tag[idx], r2rmsa->gf_tag[idx+1]);
	esl_sprintf(&r2rmsa->gf[idx],     r2rmsa->gf[idx+1]);
      }
      r2rmsa->ngf --;
    }
    
    esl_msa_AddGF(r2rmsa, "R2R keep all", -1, "", -1);
  }
  else {
    for (tagidx = 0; tagidx < r2rmsa->ngf; tagidx++) {
      esl_strchop(r2rmsa->gf[tagidx], -1);
      if (strcmp(r2rmsa->gf[tagidx], "keep allpairs") == 0) break;
    }
    if (tagidx < r2rmsa->ngf) { //remove 
      for (idx = tagidx; idx < r2rmsa->ngf-1; idx++) {
	free(r2rmsa->gf_tag[idx]); r2rmsa->gf_tag[idx] = NULL; 
	free(r2rmsa->gf[idx]);      r2rmsa->gf[idx] = NULL; 
	esl_sprintf(&r2rmsa->gf_tag[idx], r2rmsa->gf_tag[idx+1]);
	esl_sprintf(&r2rmsa->gf[idx],     r2rmsa->gf[idx+1]);
      }
      r2rmsa->ngf --;
    }
    esl_msa_AddGF(r2rmsa, "R2R keep allpairs", -1, "", -1);
  }
  
  /* replace the r2r 'cov_SS_cons' GC line(s) with our own */
  if (!do_r2rcovmarkup) {
    for (s = 0; s < nct; s ++) {
      if (s == 0) esl_sprintf(&covtag1, "%s", covtag);
      else        esl_sprintf(&covtag1, "%s_%d", covtag, s);
      for (tagidx = 0; tagidx < r2rmsa->ngc; tagidx++)
	if (strcmp(r2rmsa->gc_tag[tagidx], covtag1) == 0) break;
      if (tagidx == r2rmsa->ngc) {
	ESL_REALLOC(r2rmsa->gc_tag, (r2rmsa->ngc+1) * sizeof(char **));
	ESL_REALLOC(r2rmsa->gc,     (r2rmsa->ngc+1) * sizeof(char **));
	r2rmsa->gc[r2rmsa->ngc] = NULL;
	r2rmsa->ngc++;
      }
      if (r2rmsa->gc_tag[tagidx]) free(r2rmsa->gc_tag[tagidx]); r2rmsa->gc_tag[tagidx] = NULL;
      if ((status = esl_strdup(covtag1, -1, &(r2rmsa->gc_tag[tagidx]))) != eslOK) goto ERROR;
      free(r2rmsa->gc[tagidx]); r2rmsa->gc[tagidx] = NULL; 
      esl_sprintf(&(r2rmsa->gc[tagidx]), "%s", prv_covstr);
    }
    if (verbose) esl_msafile_Write(stdout, r2rmsa, eslMSAFILE_PFAM);
  }

  // 'SS_cons_<n>' are potential pseudoknots draw them
  // this is the R2R incantation to draw an outline around them
  r2r_pseudoknot(r2rmsa, nct, ctlist);
  
 /* write the R2R annotated to PFAM format */
  if ((fp = fopen(r2rfile, "w")) == NULL) esl_fatal("Failed to open r2rfile %s", r2rfile);
  esl_msafile_Write(fp, r2rmsa, eslMSAFILE_PFAM);
  fclose(fp);

  /* produce the R2R pdf */
  if (makepdf) {
    status = r2r_depict_pdf(r2rfile, verbose, errbuf);
    if (status != eslOK) goto ERROR;
  }
  if (makesvg) {
    status = r2r_depict_svg(r2rfile, verbose, errbuf);
    if (status != eslOK) goto ERROR;
  }

  esl_msa_Destroy(r2rmsa);

  remove(tmpinfile);
  remove(tmpoutfile);
  
  free(tok);
  if (tag) free(tag);
  if (cmd)  free(cmd);
  if (args) free(args);
  if (ssstr) free(ssstr);
  if (covstr) free(covstr);
  if (prv_covstr) free(prv_covstr);
  return eslOK;

 ERROR:
  remove(tmpinfile);
  remove(tmpoutfile);
  
  if (msa)    esl_msa_Destroy(msa);
  if (r2rmsa) esl_msa_Destroy(r2rmsa);
  if (tag)    free(tag);
  if (tok)    free(tok);
  if (cmd)    free(cmd);
  if (args)   free(args);
  if (ssstr)  free(ssstr);
  if (covstr) free(covstr);
  if (prv_covstr) free(prv_covstr);
  return status;
}


/*------- internal functions -------------------*/
static int
r2r_depict_pdf(char *r2rfile, int verbose, char *errbuf)
{
  char *r2rpdf = NULL;
  char *args = NULL;
  char *cmd = NULL;
  int   status;

  /* produce the R2R pdf */
  if (RSCAPE_BIN)                         // look for the installed executable
    esl_sprintf(&cmd, "%s/r2r", RSCAPE_BIN);  
  else
    ESL_XFAIL(status, errbuf, "Failed to find R2R executable\n");
  
  esl_sprintf(&r2rpdf, "%s.pdf", r2rfile);
  esl_sprintf(&args, "%s %s %s >/dev/null", cmd, r2rfile, r2rpdf);
 
  status = system(args);
  if (status == -1) ESL_XFAIL(status, errbuf, "Failed to run R2R2pdf\n");

  free(cmd);
  free(args);
  free(r2rpdf);
  
  return eslOK;

 ERROR:
  if (cmd)  free(cmd);
  if (args) free(args);
  return status;
 }


static int
r2r_depict_svg(char *r2rfile, int verbose, char *errbuf)
{
  char *r2rsvg = NULL;
  char *args = NULL;
  char *cmd = NULL;
  int   status;

  /* produce the R2R svg */
  if (RSCAPE_BIN)  // look for the installed executable
    esl_sprintf(&cmd, "%s/r2r", RSCAPE_BIN);  
  else
    ESL_XFAIL(status, errbuf, "Failed to find R2R executable\n");

  esl_sprintf(&r2rsvg, "%s.svg", r2rfile);
  esl_sprintf(&args, "%s %s %s >/dev/null", cmd, r2rfile, r2rsvg);
  status = system(args);
  if (status == -1) ESL_XFAIL(status, errbuf, "Failed to run R2Rsvg\n");

  if (cmd)    free(cmd);
  if (args)   free(args);
  if (r2rsvg) free(r2rsvg);
  return eslOK;
 
 ERROR:
  if (cmd)    free(cmd);
  if (args)   free(args);
  if (r2rsvg) free(r2rsvg);
  return status;
}


static int
r2r_pseudoknot(ESL_MSA *msa, int nct, int **ctlist)
{
  char       *ss   = NULL;
  char       *new  = NULL;
  char       *tag  = NULL;
  char        ssi;
  int         L = msa->alen;
  int         i;
  int         s;
  int         status;
  
  // s= 0 correspond to the main nested structure, all the other
  // we annotate a pseudoknots

  
  /* Initialization */
  ESL_ALLOC(ss,  sizeof(char) * (L+1));
  ESL_ALLOC(new, sizeof(char) * (L+1));
  
  for (s = 1; s < nct; s ++) {
    esl_ct2wuss(ctlist[s], L, ss);

    // mark the whole region to form a callout
    new[L] = '\0';
    for (i = 1; i <= L; i ++) {
      ssi = ss[i-1];
      
      // mark the region. It should be a hairpin possibly with bulges and multiloops
      // but no other structure is allowed 
      if      (ssi == '<' ||  ssi == '>' ||  ssi == '-') new[i-1] = '.'; // a helix or bulge or internal loop
      else if (ssi == '_')                               new[i-1] = '-'; // the hairpin loop
      else if (ssi == ':')                               new[i-1] = '-'; // everything else should be unpaired
      else  status = eslFAIL;
    }
    
    esl_sprintf(&tag, "SUBFAM_pknot%d_R2R_LABEL", s);
    esl_msa_AppendGC(msa, tag, new);
    
    // other markups necessary for r2r to plot the pseudokntos
    esl_sprintf(&tag, "R2R ignore_ss_except_for_pairs _%d outline", s);
    esl_msa_AddGF(msa, tag, -1, "", -1);

    esl_sprintf(&tag, "SUBFAM_PERL_PRED pknot%d return 1;", s);
    esl_msa_AddGF(msa, tag, -1, "", -1);
    esl_sprintf(&tag, "SUBFAM_pknot%d_R2R subst_ss _%d primary", s, s);
    esl_msa_AddGF(msa, tag, -1, "", -1);
    esl_sprintf(&tag, "SUBFAM_pknot%d_R2R no5", s);
    esl_msa_AddGF(msa, tag, -1, "", -1);
    esl_sprintf(&tag, "SUBFAM_pknot%d_R2R outline_nuc all", s);
    esl_msa_AddGF(msa, tag, -1, "", -1);
    esl_sprintf(&tag, "SUBFAM_pknot%d_R2R set_dir pos0 90 f", s);
    esl_msa_AddGF(msa, tag, -1, "", -1);
  }

  free(ss);
  free(new);
  free(tag);
  
  return eslOK;
  
 ERROR:
  if (ss)  free(ss);
  if (new) free(new);
  if (tag) free(tag);
  return status;
}
