/* covariation.c */

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
#include "esl_exponential.h"
#include "esl_gamma.h"
#include "esl_histogram.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "contactmap.h"
#include "covariation.h"
#include "correlators.h"
#include "cococyk.h"
#include "covgrammars.h"
#include "cykcov.h"
#include "logsum.h"
#include "pottsbuild.h"
#include "pottsscore.h"
#include "power.h"
#include "r2rdepict.h"
#include "ratematrix.h"
#include "ribosum_matrix.h"
#include "structure.h"

static double cov2evalue(double cov, int Nc, ESL_HISTOGRAM *h, double *surv);
static double evalue2cov(double eval, int Nc, ESL_HISTOGRAM *h, double *survfit);
static double cov_histogram_pmass(ESL_HISTOGRAM *h, double target_pmass, double target_fracfit);
static int    cov_histogram_plotdensity(FILE *pipe, ESL_HISTOGRAM *h, double *survfit, char *key, double posx, double posy, int logval, int style1, int style2);
static int    cov_histogram_plotsurvival(FILE *pipe, ESL_HISTOGRAM *h, double *survfit, char *key, double posx, double posy, int logval, int style1, int style2);
static int    cov_histogram_plotexpectsurv(FILE *pipe, int Nc, ESL_HISTOGRAM *h, double *survfit, char *key, char *axes, int addtics, double posx, double posy, int logval, 
					   int linespoints, int style1, int style2);
static int    cov_histogram_plotqq(FILE *pipe, struct data_s *data, ESL_HISTOGRAM *h1, ESL_HISTOGRAM *h2, char *key, int logval, 
				   int linespoints, int style1, int style2);
static int    cov_plot_lineatexpcov(FILE *pipe, struct data_s *data, double expsurv, int Nc, ESL_HISTOGRAM *h, double *survfit, ESL_HISTOGRAM *h2,

				    
				    char *axes, char *key, double ymax, double ymin, double xmax, double xmin, int style1, int style2);
static int    cov_plot_extra_yaxis(FILE *pipe, double ymax, double ymin, double xoff, char *ylabel, int style);

int                 
cov_Calculate(struct data_s *data, ESL_MSA *msa, RANKLIST **ret_ranklist, HITLIST **ret_hitlist, int analyze)
{
  RANKLIST      *ranklist = NULL;
  HITLIST       *hitlist = NULL;
  int          **ctlist = NULL;
  char          *ss = NULL;
  COVCLASS       covclass = data->mi->class;
  int            shiftnonneg = FALSE;
  int            nct = 0;
  int            s;
  int            status;

  /* Calculate the covariation matrix */
  switch(data->covmethod) {
  case NONPARAM:
  case AKMAEV:
  case POTTS:
    if ( !(data->covtype == RAF  || data->covtype == RAFp  || data->covtype == RAFa ||
	   data->covtype == RAFS || data->covtype == RAFSp || data->covtype == RAFSa ) ) {
      status = corr_Probs(data->r, msa, data->T, data->ribosum, data->mi, data->covmethod, data->tol, data->verbose, data->errbuf);
      if (status != eslOK) goto ERROR;
    }
    break;
  }
 
  switch(data->covtype) {
  case PTFp:
  case PTAp:
  case PTDp:
    shiftnonneg = TRUE; // gremling shifts scores to be nonnegative
    status = potts_CalculateCOV(data);
    if (status != eslOK) goto ERROR; 
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR;
    break;
  case CHIa: 
    status = corr_CalculateCHI         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case CHIp:
    status = corr_CalculateCHI         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR;  
    break;
  case CHI: 
    status = corr_CalculateCHI         (covclass, data, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;     
    break;
  case GTa: 
    status = corr_CalculateGT          (covclass, data, FALSE,  NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case GTp: 
    status = corr_CalculateGT          (covclass, data, FALSE,  NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case GT: 
    status = corr_CalculateGT          (covclass, data, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case MIa: 
    status = corr_CalculateMI          (covclass, data, FALSE,  NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case MIp: 
    status = corr_CalculateMI          (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case MI: 
    status = corr_CalculateMI          (covclass, data, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case MIra: 
    status = corr_CalculateMIr         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case MIrp:
    status = corr_CalculateMIr         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR;  
    break;
  case MIr: 
    status = corr_CalculateMIr         (covclass, data, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case MIga: 
    status = corr_CalculateMIg         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case MIgp:
    status = corr_CalculateMIg         (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR;  
    break;
  case MIg: 
    status = corr_CalculateMIg          (covclass, data, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case OMESa: 
    status = corr_CalculateOMES        (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case OMESp: 
    status = corr_CalculateOMES        (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case OMES: 
    status = corr_CalculateOMES        (covclass, data, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
 case RAFa: 
   status = corr_CalculateRAF         (covclass, data, msa, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case RAFp: 
    status = corr_CalculateRAF         (covclass, data, msa, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case RAF:
    status = corr_CalculateRAF         (covclass, data, msa, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case RAFSa: 
    status = corr_CalculateRAFS        (covclass, data, msa, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case RAFSp: 
    status = corr_CalculateRAFS        (covclass, data, msa, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,      data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case RAFS: 
    status = corr_CalculateRAFS        (covclass, data, msa, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  case CCFa: 
    status = corr_CalculateCCF        (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(ASC,     data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case CCFp: 
    status = corr_CalculateCCF        (covclass, data, FALSE,   NULL,      NULL);
    if (status != eslOK) goto ERROR;
    status = corr_CalculateCOVCorrected(APC,     data, analyze, &ranklist, &hitlist, shiftnonneg);
    if (status != eslOK) goto ERROR; 
    break;
  case CCF: 
    status = corr_CalculateCCF        (covclass, data, analyze, &ranklist, &hitlist);
    if (status != eslOK) goto ERROR;
    break;
  default:
    ESL_XFAIL(eslFAIL, data->errbuf, "wrong covariation type\n");
    break;
  }
  if (data->mode != RANSS) fprintf(data->sumfp, "\n");   

  if (!data->nofigures && data->mode == GIVSS && data->nbpairs > 0 && hitlist) { // do the plots only for GIVSS
    if  (msa->abc->type == eslRNA || msa->abc->type == eslDNA) {
      struct_SplitCT(data->ct, msa->alen, &nct, &ctlist, TRUE);
      
      status = struct_DotPlot(data->gnuplot, data->dplotfile, msa, nct, ctlist, data->mi, data->msamap, data->firstpos, data->samplesize, hitlist,
			   TRUE, data->verbose, data->errbuf);
      if  (status != eslOK) goto ERROR;
      status = struct_DotPlot(data->gnuplot, data->dplotfile, msa, nct, ctlist, data->mi, data->msamap, data->firstpos, data->samplesize, hitlist,
			   FALSE, data->verbose, data->errbuf);
      if  (status != eslOK) goto ERROR;
      status = r2r_Depict(data->R2Rfile, data->R2Rall, msa, nct, ctlist, hitlist, TRUE, TRUE, data->verbose, data->errbuf);
      if  (status != eslOK) goto ERROR;
    }
  }
  
  if (ret_ranklist) *ret_ranklist = ranklist; else if (ranklist) cov_FreeRankList(ranklist);
  if (ret_hitlist)  *ret_hitlist = hitlist;   else if (hitlist)  cov_FreeHitList(hitlist);
  for (s = 0; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  if (ctlist) free(ctlist);
  return eslOK;
  
 ERROR:
  if (ranklist) cov_FreeRankList(ranklist);
  if (hitlist)  cov_FreeHitList(hitlist);
  for (s = 0; s < nct; s ++) if (ctlist[s]) free(ctlist[s]);
  if (ctlist) free(ctlist);
  return status;
}


int 
cov_THRESHTYPEString(char **ret_threshtype, THRESHTYPE type, char *errbuf)
{
  int status;

  switch(type) {
  case Eval:     esl_sprintf(ret_threshtype, "Eval");    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "wrong THRESHTYPE");
  }

  return eslOK;
  
 ERROR:
  return status;
}

int
cov_SignificantPairs_Ranking(struct data_s *data, RANKLIST **ret_ranklist, HITLIST **ret_hitlist)
{
  struct mutual_s *mi  = data->mi;
  ESL_DMATRIX     *mtx = mi->COV;
  char            *threshtype = NULL;
  char            *covtype = NULL;
  RANKLIST        *ranklist = NULL;
  HITLIST         *hitlist = NULL;
  double           pmass, newmass;
  double           bmax;
  double           add;
  double           cov;
  double           eval;
  double           tol = 1e-2;
  int              select;
  int              i, j;
  int              b;
  int              status;

  corr_COVTYPEString(&covtype, mi->type, data->errbuf);
  cov_THRESHTYPEString(&threshtype, data->thresh->type, NULL);

  // really bad alignment, covariation spread is smaller than tol, stop here.
  if (data->w < tol) {
    if (data->mode == GIVSS) {
      fprintf(stdout,    "#\n# Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] \n");
      fprintf(stdout,    "# %s    %g           [%.2f,%.2f]    [%d | %d %d %d | %.2f %.2f %.2f] \n#\n", 
	      covtype, data->thresh->val, mi->minCOV, mi->maxCOV, 0, 0, data->clist->ncnt, 0, 0.0, 0.0, 0.0);    
      fprintf(stdout, "#-------------------------------------------------------------------------------------------------------\n");
      fprintf(stdout, "covariation scores are almost constant, no further analysis.\n");
    }

    free(threshtype); 
    free(covtype);
    return eslOK; 
  }
  
  /* histogram parameters */
  bmax = mi->maxCOV + 5*data->w;
  while (fabs(bmax-data->bmin) < tol) bmax += data->w;
  ranklist = cov_CreateRankList(bmax, data->bmin, data->w);
 
  // Create the histograms
  for (i = 0; i < mi->alen-1; i ++) 
    for (j = i+1; j < mi->alen; j ++) {

      // if found in the pdb structure, make sure
      // they are at least mind appart
      if (data->msa2pdb[i] >= 0 &&
	  data->msa2pdb[j] >= 0 &&
	  data->msa2pdb[j]-data->msa2pdb[i] < data->clist->mind)
	continue;
      
      /* add to the ha histogram  */
      add = ESL_MAX(mtx->mx[i][j], data->bmin+data->w);
      esl_histogram_Add(ranklist->ha, add);
 
      /* add to the histogram of base pairs (hb) and not bps (ht) */
      if (data->mode == GIVSS || data->mode == CYKSS) {

	switch(data->samplesize) {
	case SAMPLE_CONTACTS:
	  select = CMAP_IsContactLocal(i+1,j+1,data->clist);
	  break;
	case  SAMPLE_BP:
	  select = CMAP_IsBPLocal(i+1,j+1,data->clist);
	  break;
	case  SAMPLE_WC:
	  select = CMAP_IsWCLocal(i+1,j+1,data->clist);
	  break;	  
	case  SAMPLE_ALL:
	  select = FALSE;
	  break;	  
	}

	if (select)
	  esl_histogram_Add(ranklist->hb, add);
	else 
	  esl_histogram_Add(ranklist->ht, add);
      }
    }

  /* Histogram and Fit */
  if (data->ranklist_null && data->mode == GIVSS ) {
    if (data->ranklist_null->ha->nb < ranklist->ha->nb) {
      ESL_REALLOC(data->ranklist_null->ha->obs, sizeof(uint64_t) * ranklist->ha->nb);
      for (i = data->ranklist_null->ha->nb; i < ranklist->ha->nb; i++) data->ranklist_null->ha->obs[i] = 0;
      data->ranklist_null->ha->nb = ranklist->ha->nb;
    }
    
    /* censor the histogram and do an exponential fit to the tail */
    pmass = cov_histogram_pmass(data->ranklist_null->ha, data->pmass, data->fracfit);
    if (data->doexpfit) {
      status = cov_NullFitExponential(data->ranklist_null->ha, &data->ranklist_null->survfit, pmass,
				      &newmass, &data->mu, &data->lambda, data->verbose, data->errbuf);
      if (status != eslOK) ESL_XFAIL(eslFAIL, data->errbuf, "bad exponential fit.");
      if (data->verbose) {
	fprintf(data->outfp, "# ExpFIT: pmass %f mu %f lambda %f\n", newmass, data->mu, data->lambda);
	fprintf(stdout,      "# ExpFIT: pmass %f mu %f lambda %f\n", newmass, data->mu, data->lambda);
      }
    }
    else { // a gamma fit
      status = cov_NullFitGamma(data->ranklist_null->ha, &data->ranklist_null->survfit, pmass,
				&newmass, &data->mu, &data->lambda, &data->tau, data->verbose, data->errbuf);      
      if (status != eslOK) ESL_XFAIL(eslFAIL, data->errbuf, "bad Gamma fit.");
      if (data->verbose) {
	fprintf(data->outfp, "# GammaFIT: pmass %f mu %f lambda %f tau %f\n", newmass, data->mu, data->lambda, data->tau);
	fprintf(stdout,      "# GammaFIT: pmass %f mu %f lambda %f tau %f\n", newmass, data->mu, data->lambda, data->tau);
      }
    }
  }
  
  // assign the covthresh corresponding to the evalue threshold
  if (data->mode == GIVSS) {
    data->thresh->sc_bp =
      (data->ranklist_null)? evalue2cov(data->thresh->val, (ranklist->hb->Nc>0)? ranklist->hb->Nc:ranklist->ha->Nc, data->ranklist_null->ha, data->ranklist_null->survfit) : eslINFINITY;
    data->thresh->sc_nbp =
      (data->ranklist_null)? evalue2cov(data->thresh->val, (ranklist->hb->Nc>0)? ranklist->ht->Nc:ranklist->ha->Nc, data->ranklist_null->ha, data->ranklist_null->survfit) : eslINFINITY;
  }

  status = cov_ROC(data, covtype, ranklist);
  if (status != eslOK) ESL_XFAIL(eslFAIL, data->errbuf, "bad ROCfile.");

  if (data->mode == GIVSS) {
    status = cov_CreateHitList(data, mi, ranklist, &hitlist, covtype, threshtype);
    if (status != eslOK) goto ERROR;
  } 
  if (ret_ranklist) *ret_ranklist = ranklist; else if (ranklist) cov_FreeRankList(ranklist);
  if (ret_hitlist)  *ret_hitlist  = hitlist;  else if (hitlist)  cov_FreeHitList(hitlist);

  if (threshtype) free(threshtype); 
  if (covtype)    free(covtype);
  
  return eslOK;
  
 ERROR:
  if (ranklist)   cov_FreeRankList(ranklist);
  if (hitlist)    cov_FreeHitList(hitlist);
  if (threshtype) free(threshtype); 
  if (covtype)    free(covtype);
  if (data->ranklist_null && (data->mode == GIVSS || data->mode == CYKSS)) free(data->ranklist_null->survfit);
  return status;
}

int
cov_ROC(struct data_s *data, char *covtype, RANKLIST *ranklist)
{
  struct mutual_s *mi = data->mi;
  ESL_DMATRIX     *mtx = mi->COV;
  ESL_DMATRIX     *eval = NULL;
  double           sen;
  double           ppv;
  double           F;
  double           target_cov;
  double           target_eval;
  double           cov;
  double           E;
  int              L = mi->alen;
  int              select;
  int              i, j;
  int              b;
  int              fp, tf, t, f, neg;

  if (data->rocfp == NULL) return eslOK;
  if (data->ranklist_null == NULL) return eslOK; 

  eval = esl_dmatrix_Create(L, L);
  for (i = 0; i < L-1; i ++) 
    for (j = i+1; j < L; j ++) {
      cov = mtx->mx[i][j];
      
      switch(data->samplesize) {
      case SAMPLE_CONTACTS:
	select =  CMAP_IsContactLocal(i+1, j+1, data->clist);
	break;
      case SAMPLE_BP:
	select = CMAP_IsBPLocal(i+1, j+1, data->clist);
	break;
      case SAMPLE_WC:
	select = CMAP_IsWCLocal(i+1, j+1, data->clist);
	break;
      case SAMPLE_ALL:
	select = FALSE;
	break;	
      }
      
      if (select) eval->mx[i][j] = cov2evalue(cov, ranklist->hb->Nc, data->ranklist_null->ha, data->ranklist_null->survfit);
      else        eval->mx[i][j] = cov2evalue(cov, ranklist->ht->Nc, data->ranklist_null->ha, data->ranklist_null->survfit);
    }    
  
  if (data->mode == GIVSS || data->mode == CYKSS) {
    fprintf(data->rocfp, "# Method: %s\n", covtype);      
    if (mi->ishuffled)
      fprintf(data->rocfp, "#shuffled_cov_score    FP              TP           Found            True       Negatives        Sen     PPV     F       E-value\n"); 
    else
      fprintf(data->rocfp, "#cov_score             FP              TP           Found            True       Negatives        Sen     PPV     F       E-value\n"); 
  }
 
  for (b = ranklist->ha->imax; b >= ranklist->ha->imin; b --) {
    target_cov  = esl_histogram_Bin2LBound(ranklist->ha, b);
    target_eval = cov2evalue(target_cov, ranklist->ha->Nc, data->ranklist_null->ha, data->ranklist_null->survfit);
    
    f = t = tf = 0;
    for (i = 0; i < L-1; i ++) 
      for (j = i+1; j < L; j ++) {
	
	// if found in the pdb structure, make sure
	// they are at least mind appart
	if (data->msa2pdb[i] >= 0 &&
	    data->msa2pdb[j] >= 0 &&
	    data->msa2pdb[j]-data->msa2pdb[i] < data->clist->mind) continue;
	
	cov = mtx->mx[i][j];
	E   = eval->mx[i][j];
	
	switch(data->samplesize) {
	case SAMPLE_CONTACTS:
	  select =  CMAP_IsContactLocal(i+1, j+1, data->clist);
	  break;
	case SAMPLE_BP:
	  select = CMAP_IsBPLocal(i+1, j+1, data->clist);
	  break;
	case SAMPLE_WC:
	  select = CMAP_IsWCLocal(i+1, j+1, data->clist);
	  break;
	case SAMPLE_ALL:
	  select = CMAP_IsContactLocal(i+1, j+1, data->clist);
	  break;	
	}
	
	if (E <= target_eval)    f  ++;
	if (select)            { t  ++;
	  if (E <= target_eval)  tf ++;
	}
     }    
    fp  = f - tf;
    sen = (t > 0)? 100. * (double)tf / (double)t : 0.0;
    ppv = (f > 0)? 100. * (double)tf / (double)f : 0.0;
    F   = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;   
    neg  = L * (L-1) / 2 - t;
    
    if (data->mode == GIVSS || data->mode == CYKSS) 
      fprintf(data->rocfp, "%.5f\t%8d\t%8d\t%8d\t%8d\t%8d\t%.2f\t%.2f\t%.2f\t%g\n", target_cov, fp, tf, f, t, neg, sen, ppv, F, target_eval);
  }

  esl_dmatrix_Destroy(eval);
  return eslOK;
}

RANKLIST *
cov_CreateRankList(double bmax, double bmin, double w)
{
  RANKLIST *ranklist = NULL;
  int       status;
  
  ESL_ALLOC(ranklist, sizeof(RANKLIST));

  ranklist->ha = NULL;
  ranklist->ht = NULL;
  ranklist->hb = NULL;
  
  ranklist->ha = esl_histogram_CreateFull(bmin, bmax, w);
  ranklist->ht = esl_histogram_CreateFull(bmin, bmax, w);
  ranklist->hb = esl_histogram_CreateFull(bmin, bmax, w);
  if (ranklist->ha == NULL) goto ERROR;
  if (ranklist->ht == NULL) goto ERROR;
  if (ranklist->hb == NULL) goto ERROR;

  ranklist->survfit = NULL;
  
  return ranklist;

 ERROR:
  return NULL;
}

int 
cov_Add2SubsHistogram(ESL_HISTOGRAM *hsubs, HITLIST *hitlist, int verbose)
{
  int h;

  for (h = 0; h < hitlist->nhit; h ++) 
    if (hitlist->hit[h].bptype == WWc) 
      esl_histogram_Add(hsubs, (double)(hitlist->hit[h].nsubs+1));
   
  if (verbose) 
    esl_histogram_Write(stdout, hsubs);
  
  return eslOK;
}

int
cov_GrowRankList(RANKLIST **oranklist, double bmax, double bmin)
{
  RANKLIST *ranklist = *oranklist;
  RANKLIST *new = NULL;
  double    new_bmin;
  int       b, newb;
  int       status;

  /* bmin has to be a w-multiple of ranklist->bin */
  new_bmin = ranklist->ha->bmin;
  if (bmin < ranklist->ha->bmin) new_bmin -= fabs(bmin) * 2. * ranklist->ha->w;
  
  new = cov_CreateRankList(ESL_MAX(bmax, ranklist->ha->bmax), new_bmin, ranklist->ha->w);
  if (new == NULL) { status = eslFAIL; goto ERROR; }

  new->ha->n    = ranklist->ha->n;
  new->ha->xmin = ranklist->ha->xmin;
  new->ha->xmax = ranklist->ha->xmax;
  new->ha->imin = ranklist->ha->imin;
  new->ha->imax = ranklist->ha->imax;
  new->ha->Nc   = ranklist->ha->Nc;
  new->ha->No   = ranklist->ha->No;

  new->ht->n    = ranklist->ht->n;
  new->ht->xmin = ranklist->ht->xmin;
  new->ht->xmax = ranklist->ht->xmax;
  new->ht->imin = ranklist->ht->imin;
  new->ht->imax = ranklist->ht->imax;
  new->ht->Nc   = ranklist->ht->Nc;
  new->ht->No   = ranklist->ht->No;

  new->hb->n    = ranklist->hb->n;
  new->hb->xmin = ranklist->hb->xmin;
  new->hb->xmax = ranklist->hb->xmax;
  new->hb->imin = ranklist->hb->imin;
  new->hb->imax = ranklist->hb->imax;
  new->hb->Nc   = ranklist->hb->Nc;
  new->hb->No   = ranklist->hb->No;

  for (b = ranklist->ha->imin; b <= ranklist->ha->imax; b ++) {
    cov_ranklist_Bin2Bin(b, ranklist->ha, new->ha, &newb);
    if (newb < new->ha->nb) {
      new->ha->obs[newb] = ranklist->ha->obs[b];
    }
  }
    
  cov_FreeRankList(ranklist);
  *oranklist = new;
  return eslOK;

 ERROR:
  if (new) cov_FreeRankList(new);
  return status;
}

int              
cov_DumpRankList(FILE *fp, RANKLIST *ranklist)
{
  double cov;
  double eval;
  int    b;

  if (ranklist == NULL) return eslOK;
  
  printf("imin %d imax %d covmin %f covmax %f\n", ranklist->ha->imin, ranklist->ha->imax, ranklist->ha->xmin, ranklist->ha->xmax);
  for (b = ranklist->ha->imax; b >= ranklist->ha->imin; b --) {
    cov = esl_histogram_Bin2LBound(ranklist->ha,b);
    eval = cov2evalue(cov, ranklist->ha->Nc, ranklist->ha, NULL);
    printf("cov %f eval %g\n", cov, eval);
  }
 
  return eslOK;
}
int              
cov_DumpHistogram(FILE *fp, ESL_HISTOGRAM *h)
{
  double nobs = 0.;
  double cum_obs = 0.;
  double cum_exp = 0.;
  int    b;

  for (b = h->imax; b >= h->imin; b --) nobs += (double)h->obs[b];

  printf("imin %d imax %d covmin %f covmax %f nobs %f\n", h->imin, h->imax, h->xmin, h->xmax, nobs);
  for (b = h->imax; b >= h->imin; b --) {
    cum_obs += (double)h->obs[b];
    if (h->expect) cum_exp += (double)h->expect[b];
    if (h->expect) printf("b %d val %f obs %f cum_obs %f expect %f cum_exp %f\n", 
			  b, esl_histogram_Bin2LBound(h, b), (double)h->obs[b], cum_obs, h->expect[b], cum_exp); 
    else           printf("b %d val %f obs %f cum_obs %f\n", b, esl_histogram_Bin2LBound(h, b), (double)h->obs[b], cum_obs); 
  }
  
  return eslOK;
}

int 
cov_CreateHitList(struct data_s *data, struct mutual_s *mi, RANKLIST *ranklist, HITLIST **ret_hitlist, char *covtype, char *threshtype)
{
  HITLIST  *hitlist = NULL;
  int      *msa2pdb = data->msa2pdb;
  double    sen, ppv, F;
  double    cov;
  double    eval;
  BPTYPE    bptype;
  int       select;
  int       is_compatible;
  int       alloc_nhit = 5;
  int       tf = 0;
  int       f  = 0;
  int       t, fp;
  int       nhit;
  int       h = 0;
  int64_t   n = 0;
  int64_t   i, j;
  int       status;
  
  ESL_ALLOC(hitlist, sizeof(HITLIST));
  hitlist->hit    = NULL;
  hitlist->srthit = NULL;

  nhit = alloc_nhit;
  ESL_ALLOC(hitlist->hit,    sizeof(HIT)   * nhit);
  ESL_ALLOC(hitlist->srthit, sizeof(HIT *) * nhit);
  hitlist->srthit[0] = hitlist->hit;

  for (i = 0; i < mi->alen-1; i++) 
    for (j = i+1; j < mi->alen; j++) {
      
      // if found in the pdb structure, make sure
      // they are at least mind appart
      if (msa2pdb[i] >= 0 &&
	  msa2pdb[j] >= 0 &&
	  msa2pdb[j]-msa2pdb[i] < data->clist->mind) continue;

      is_compatible = FALSE;
      bptype        = CMAP_GetBPTYPE(i+1, j+1, data->clist);
     
      if (data->abcisRNA && data->ct[i+1] == 0 && data->ct[j+1] == 0) {
	is_compatible = TRUE;
      }

      switch(data->samplesize) {
      case SAMPLE_CONTACTS:
	select = (bptype < BPNONE)?  TRUE : FALSE;
	break;
      case SAMPLE_BP:
	select = (bptype < STACKED)? TRUE : FALSE;
	break;
      case SAMPLE_WC:
	select = (bptype == WWc)?    TRUE : FALSE;
	break;
      case SAMPLE_ALL:
	select = FALSE;
	break;	
      }
      
      cov = mi->COV->mx[i][j];
      if (data->ranklist_null == NULL) eval = 0; // naive method: eval does not inform of statistical significance
      else {
	if (select)
	  eval = cov2evalue(cov, ranklist->hb->Nc, data->ranklist_null->ha, data->ranklist_null->survfit);
	else 
	  eval = cov2evalue(cov, ranklist->ht->Nc, data->ranklist_null->ha, data->ranklist_null->survfit);
      }

      if (eval < data->thresh->val) {	
	if (h == nhit - 1) {
 	  nhit += alloc_nhit;
	  
	  ESL_REALLOC(hitlist->hit,    sizeof(HIT)   * nhit);
	  ESL_REALLOC(hitlist->srthit, sizeof(HIT *) * nhit);
	}

	/* assign */
	hitlist->hit[h].i             = i;
	hitlist->hit[h].j             = j;
	hitlist->hit[h].sc            = cov;
	hitlist->hit[h].Eval          = eval;
	hitlist->hit[h].nsubs         = (data->spair)? data->spair[n].nsubs : 0;
	hitlist->hit[h].power         = (data->spair)? data->spair[n].power : 0.;
	hitlist->hit[h].bptype        = bptype;
	hitlist->hit[h].is_compatible = is_compatible;
	h ++;
      }
 
      n ++;
    }
  nhit = h;
  hitlist->nhit = nhit;

  switch(data->samplesize) {
  case SAMPLE_CONTACTS:
    t = data->clist->ncnt;
    break;
  case SAMPLE_BP:
    t = data->clist->nbps;
    break;
  case SAMPLE_WC:
    t = data->clist->nwwc;
    break;
  case SAMPLE_ALL: // use all contacts here
    t = data->clist->ncnt;
    break;	
  }
  
  for (h = 0; h < nhit; h ++) {
    f ++;
    
    switch(data->samplesize) {
    case SAMPLE_CONTACTS:
      select = (hitlist->hit[h].bptype < BPNONE)?  TRUE : FALSE;
      break;
    case SAMPLE_BP:
      select = (hitlist->hit[h].bptype < STACKED)? TRUE : FALSE;
      break;
    case SAMPLE_WC:
     select = (hitlist->hit[h].bptype == WWc)?     TRUE : FALSE;
      break;
    case SAMPLE_ALL: // use all contacts here
      select = (hitlist->hit[h].bptype < BPNONE)?  TRUE : FALSE;
     break;	
    }
    if (select) tf ++;
  }
  
  fp = f - tf;
  sen = (t > 0)? 100. * (double)tf / (double)t : 0.0;
  ppv = (f > 0)? 100. * (double)tf / (double)f : 0.0;
  F   = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;   
  
  if (data->mode != RANSS && data->sumfp) {
    fprintf(data->sumfp, " %s %d %d %d %.2f %.2f ", 
	    covtype, tf, t, f, (t > 0)? 100.*(double)tf/(double)t:0.0, (f>0)? 100.*(double)tf/(double)f:0.0);
  }

  if (data->spair != NULL && data->samplesize != SAMPLE_ALL) printf("# BPAIRS observed to covary %d\n#\n", tf);
  
  if (data->outfp) {
    CMAP_DumpShort(data->outfp, data->clist);
    fprintf(data->outfp,    "#\n# Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] \n");
    fprintf(data->outfp,    "# %s    %g           [%.2f,%.2f]    [%d | %d %d %d | %.2f %.2f %.2f] \n#\n", 
	    covtype, data->thresh->val, ranklist->ha->xmin, ranklist->ha->xmax, fp, tf, t, f, sen, ppv, F);
    cov_WriteHitList(data->outfp, nhit, hitlist, data->msamap, data->firstpos);
  }
  
  if (data->outsrtfp) {
    CMAP_DumpShort(data->outsrtfp, data->clist);
    fprintf(stdout,         "#\n# Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] \n");
    fprintf(data->outsrtfp, "#\n# Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] \n");
    fprintf(stdout,         "# %s    %g         [%.2f,%.2f]     [%d | %d %d %d | %.2f %.2f %.2f] \n#\n", 
	    covtype, data->thresh->val, ranklist->ha->xmin, ranklist->ha->xmax, fp, tf, t, f, sen, ppv, F);
    fprintf(data->outsrtfp, "# %s    %g         [%.2f,%.2f]     [%d | %d %d %d | %.2f %.2f %.2f] \n#\n", 
	    covtype, data->thresh->val, ranklist->ha->xmin, ranklist->ha->xmax, fp, tf, t, f, sen, ppv, F);
    cov_WriteRankedHitList(stdout,         nhit, hitlist, data->msamap, data->firstpos, data->statsmethod);
    cov_WriteRankedHitList(data->outsrtfp, nhit, hitlist, data->msamap, data->firstpos, data->statsmethod);
  }
  
  if (ret_hitlist) *ret_hitlist = hitlist; else cov_FreeHitList(hitlist);
  return eslOK;
  
 ERROR:
  if (hitlist) cov_FreeHitList(hitlist);
  return status;
}


int 
cov_CreateCYKHitList(struct data_s *data, RANKLIST *ranklist, HITLIST *hitlist, HITLIST **ret_cykhitlist, char *covtype, char *threshtype)
{
  HITLIST  *cykhitlist = NULL;
  SPAIR    *spair      = NULL;
  double    sen, ppv, F;
  int       dim;
  int       nhit = (hitlist)? hitlist->nhit : 0;
  int       select;
  int       tf = 0;
  int       f  = 0;
  int       t, fp;
  int       h = 0;
  int       status;

  if (nhit > 0) {
    ESL_ALLOC(cykhitlist, sizeof(HITLIST));
    cykhitlist->hit    = hitlist->hit;
    cykhitlist->srthit = hitlist->srthit;
    
    cykhitlist->nhit = nhit;
    ESL_ALLOC(cykhitlist->hit,    sizeof(HIT)   * nhit);
    ESL_ALLOC(cykhitlist->srthit, sizeof(HIT *) * nhit);
    cykhitlist->srthit[0] = cykhitlist->hit;
    
    // the cyk hitlist is identical expect for the basepair annotation
    for (h = 0; h < nhit; h ++) {
      cykhitlist->hit[h].i             = hitlist->hit[h].i;
      cykhitlist->hit[h].j             = hitlist->hit[h].j;
      cykhitlist->hit[h].sc            = hitlist->hit[h].sc;
      cykhitlist->hit[h].Eval          = hitlist->hit[h].Eval;
      cykhitlist->hit[h].nsubs         = hitlist->hit[h].nsubs;
      cykhitlist->hit[h].power         = hitlist->hit[h].power;
      cykhitlist->hit[h].bptype        = CMAP_GetBPTYPE(hitlist->hit[h].i+1, hitlist->hit[h].j+1, data->clist);
      cykhitlist->hit[h].is_compatible = FALSE;
      if (data->abcisRNA && data->ct[hitlist->hit[h].i+1] == 0 && data->ct[hitlist->hit[h].j+1] == 0) {
	cykhitlist->hit[h].is_compatible  = TRUE;
      }
    }
  }
  
  switch(data->samplesize) {
  case SAMPLE_CONTACTS:
    t = data->clist->ncnt;
    break;
  case SAMPLE_BP:
    t = data->clist->nbps;
    break;
  case SAMPLE_WC:
    t = data->clist->nwwc;
    break;
  case SAMPLE_ALL: // use all contacts here
    t = data->clist->ncnt;
    break;	
  }
  
  for (h = 0; h < nhit; h ++) {
    f ++;
    
    switch(data->samplesize) {
    case SAMPLE_CONTACTS:
      select = (cykhitlist->hit[h].bptype < BPNONE)?  TRUE : FALSE;
      break;
    case SAMPLE_BP:
      select = (cykhitlist->hit[h].bptype < STACKED)? TRUE : FALSE;
      break;
    case SAMPLE_WC:
     select = (cykhitlist->hit[h].bptype == WWc)?     TRUE : FALSE;
      break;
    case SAMPLE_ALL: // use all contacts here
      select = (cykhitlist->hit[h].bptype < BPNONE)?  TRUE : FALSE;
     break;	
    }
    if (select) tf ++;
  }
  fp = f - tf;
  sen = (t > 0)? 100. * (double)tf / (double)t : 0.0;
  ppv = (f > 0)? 100. * (double)tf / (double)f : 0.0;
  F   = (sen+ppv > 0.)? 2.0 * sen * ppv / (sen+ppv) : 0.0;   
  
  if (data->mode != RANSS && data->sumfp) {
    fprintf(data->sumfp, " %s %d %d %d %.2f %.2f ", 
	    covtype, tf, t, f, (t > 0)? 100.*(double)tf/(double)t:0.0, (f>0)? 100.*(double)tf/(double)f:0.0);
  }
  
  status = power_SPAIR_Create(&dim, &spair, data->mi->alen, data->msamap, data->power, data->clist, data->nsubs, data->ndouble, data->errbuf, data->verbose);
  if (status != eslOK) ESL_XFAIL(status, data->errbuf, "%s\n", data->errbuf);
  
  if (data->outfp) {
    fprintf(data->outfp, "\n# The predicted cyk-cov structure\n");

    power_SPAIR_Write(data->outfp, dim, spair);
    
    fprintf(data->outfp,    "#\n# Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] \n");
    fprintf(data->outfp,    "# %s    %g           [%.2f,%.2f]    [%d | %d %d %d | %.2f %.2f %.2f] \n#\n", 
	    covtype, data->thresh->val, (ranklist)?ranklist->ha->xmin:0, (ranklist)?ranklist->ha->xmax:0, fp, tf, t, f, sen, ppv, F);
    cov_WriteCYKHitList(data->outfp, nhit, hitlist, cykhitlist, data->msamap, data->firstpos);
  }
  
  if (data->outsrtfp) {
    fprintf(stdout,         "\n# The predicted cyk-cov structure\n");
    fprintf(data->outsrtfp, "\n# The predicted cyk-cov structure\n");
    
    power_SPAIR_Write(stdout,         dim, spair);
    power_SPAIR_Write(data->outsrtfp, dim, spair);
   
    fprintf(stdout,         "# BPAIRS observed to covary %d\n#\n", tf);
    fprintf(data->outsrtfp, "# BPAIRS observed to covary %d\n#\n", tf);
 
    fprintf(stdout,         "#\n# Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] \n");
    fprintf(data->outsrtfp, "#\n# Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] \n");
    fprintf(stdout,         "# %s    %g         [%.2f,%.2f]     [%d | %d %d %d | %.2f %.2f %.2f] \n#\n", 
	    covtype, data->thresh->val, (ranklist)?ranklist->ha->xmin:0, (ranklist)?ranklist->ha->xmax:0, fp, tf, t, f, sen, ppv, F);
    fprintf(data->outsrtfp, "# %s    %g         [%.2f,%.2f]     [%d | %d %d %d | %.2f %.2f %.2f] \n#\n", 
	    covtype, data->thresh->val, (ranklist)?ranklist->ha->xmin:0, (ranklist)?ranklist->ha->xmax:0, fp, tf, t, f, sen, ppv, F);

    cov_WriteCYKRankedHitList(stdout,         nhit, hitlist, cykhitlist, data->msamap, data->firstpos, data->statsmethod);
    cov_WriteCYKRankedHitList(data->outsrtfp, nhit, hitlist, cykhitlist, data->msamap, data->firstpos, data->statsmethod);
  }

  if (spair) free(spair);
  if (ret_cykhitlist) *ret_cykhitlist = cykhitlist; else cov_FreeHitList(cykhitlist);
  return eslOK;
  
 ERROR:
  if (spair) free(spair);
  if (cykhitlist) cov_FreeHitList(cykhitlist);
  return status;
}


int
cov_WriteHitList(FILE *fp, int nhit, HITLIST *hitlist, int *msamap, int firstpos)
{
  int h;
  int ih, jh;

  if (fp == NULL) return eslOK;

  if (nhit == 0) {
    fprintf(fp, "#-------------------------------------------------------------------------------------------------------\n");
    fprintf(fp, "no significant pairs\n");
  }
  else {
    fprintf(fp, "#       left_pos       right_pos        score          E-value       substitutions      power\n");
    fprintf(fp, "#-------------------------------------------------------------------------------------------------------\n");
  }

  for (h = 0; h < nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    
    if (hitlist->hit[h].bptype == WWc)      { 
      fprintf(fp, "*\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power); 
    }
    else if (hitlist->hit[h].bptype < STACKED) { 
      fprintf(fp, "**\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power); 
    }
    else if (hitlist->hit[h].bptype < BPNONE && hitlist->hit[h].is_compatible) { 
      fprintf(fp, "c~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power); 
    }
   else if (hitlist->hit[h].bptype < BPNONE) { 
      fprintf(fp, "c\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power); 
    }
    else if (hitlist->hit[h].is_compatible) { 
      fprintf(fp, "~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power); 
    }
    else { 
      fprintf(fp, " \t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n",
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power); 
    }  
  }

  return eslOK;
}

int 
cov_WriteCYKHitList(FILE *fp, int nhit, HITLIST *hitlist, HITLIST *cykhitlist, int *msamap, int firstpos)
{
  int h;
  int ih, jh;

  if (!fp)         return eslOK;
  if (!cykhitlist) return eslOK;

  fprintf(fp, "# in_cyk  in_given   left_pos       right_pos      score           E-value substitutions    power\n");
  fprintf(fp, "#-------------------------------------------------------------------------------------------------------\n");
   if (nhit == 0) fprintf(fp, "no significant pairs\n");

  for (h = 0; h < nhit; h ++) {
    ih = hitlist->hit[h].i;
    jh = hitlist->hit[h].j;
    
    if (hitlist->hit[h].bptype == WWc)      {
      if (cykhitlist->hit[h].bptype == WWc) 
	fprintf(fp, "*\t*\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else if (cykhitlist->hit[h].is_compatible)
	fprintf(fp, "~\t*\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else 
	fprintf(fp, " \t*\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      }
    
    else if (hitlist->hit[h].bptype < STACKED) {
      if (cykhitlist->hit[h].bptype == WWc) 
	fprintf(fp, "*\t**\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else if (cykhitlist->hit[h].is_compatible)
	fprintf(fp, "~\t**\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else 
	fprintf(fp, " \t**\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
    }
    
    else if (hitlist->hit[h].bptype < BPNONE && hitlist->hit[h].is_compatible) {
      if (cykhitlist->hit[h].bptype == WWc) 
	fprintf(fp, "*\tc~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else if (cykhitlist->hit[h].is_compatible)
	fprintf(fp, "~\tc~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else 
	fprintf(fp, " \tc~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
    }
    
   else if (hitlist->hit[h].bptype < BPNONE) {
           if (cykhitlist->hit[h].bptype == WWc) 
	fprintf(fp, "*\tc\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else if (cykhitlist->hit[h].is_compatible)
	fprintf(fp, "~\tc\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else 
	fprintf(fp, " \tc\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
   }
    
    else if (hitlist->hit[h].is_compatible) {
            if (cykhitlist->hit[h].bptype == WWc) 
	fprintf(fp, "*\t~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else if (cykhitlist->hit[h].is_compatible)
	fprintf(fp, "~\t~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else 
	fprintf(fp, " \t~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
    }
    
    else {
            if (cykhitlist->hit[h].bptype == WWc) 
	fprintf(fp, "*\t \t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else if (cykhitlist->hit[h].is_compatible)
	fprintf(fp, "~\t \t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
      else 
	fprintf(fp, " \t \t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->hit[h].sc, hitlist->hit[h].Eval, hitlist->hit[h].nsubs, hitlist->hit[h].power);
    }  
  }

  return eslOK;
}

static int
hit_sorted_by_score(const void *vh1, const void *vh2)
{
  HIT *h1 = *((HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
  HIT *h2 = *((HIT **) vh2);

  if      (h1->sc < h2->sc) return  1;
  else if (h1->sc > h2->sc) return -1;
  else {
 
    /* report first pair first */
    int dir1 = (h1->i < h2->i ? 1 : -1);
    int dir2 = (h1->j < h2->j ? 1 : -1);
    if (dir1 != dir2) return dir2; // so if dir1 is pos (1), and dir2 is neg (-1), this will return -1, placing h1 before h2;  otherwise, vice versa
    else              return dir1;

  }
}

static int
hit_sorted_by_eval(const void *vh1, const void *vh2)
{
  HIT *h1 = *((HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
  HIT *h2 = *((HIT **) vh2);

  if      (h1->Eval > h2->Eval) return  1;
  else if (h1->Eval < h2->Eval) return -1;
  else {
 
    /* report first pair first */
    int dir1 = (h1->i < h2->i ? 1 : -1);
    int dir2 = (h1->j < h2->j ? 1 : -1);
    if (dir1 != dir2) return dir2; // so if dir1 is pos (1), and dir2 is neg (-1), this will return -1, placing h1 before h2;  otherwise, vice versa
    else              return dir1;

  }
}

int 
cov_WriteRankedHitList(FILE *fp, int nhit, HITLIST *hitlist, int *msamap, int firstpos, STATSMETHOD statsmethod)
{
  int h;
  int ih, jh;

  if (fp == NULL) return eslOK;

  for (h = 0; h < nhit; h++) hitlist->srthit[h] = hitlist->hit + h;
  if (nhit > 1) qsort(hitlist->srthit, nhit, sizeof(HIT *), (statsmethod == NAIVE)? hit_sorted_by_score:hit_sorted_by_eval);

  fprintf(fp, "#       left_pos       right_pos        score          E-value       substitutions      power\n");
  fprintf(fp, "#-------------------------------------------------------------------------------------------------------\n");
  if (nhit == 0) fprintf(fp, "no significant pairs\n");

  for (h = 0; h < nhit; h ++) {
    ih = hitlist->srthit[h]->i;
    jh = hitlist->srthit[h]->j;
    
    if (hitlist->srthit[h]->bptype == WWc) { 
      fprintf(fp, "*\t%8d\t%8d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power); 
    }
    else if (hitlist->srthit[h]->bptype < STACKED) { 
      fprintf(fp, "**\t%8d\t%8d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power); 
    }
    else if (hitlist->srthit[h]->bptype < BPNONE && hitlist->srthit[h]->is_compatible) { 
      fprintf(fp, "c~\t%8d\t%8d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power); 
    }
    else if (hitlist->srthit[h]->bptype < BPNONE) { 
      fprintf(fp, "c\t%8d\t%8d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power); 
    }
    else if (hitlist->srthit[h]->is_compatible) { 
      fprintf(fp, "~\t%8d\t%8d\t%.5f\t%g\t%lld\t\t%.2f\n", 
	      msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power); 
    }
    else { 
      fprintf(fp, " \t%8d\t%8d\t%.5f\t%g\t%lld\t\t%.2f\n",
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power); 
    }  
  }

  return eslOK;
}

int 
cov_WriteCYKRankedHitList(FILE *fp, int nhit, HITLIST *hitlist, HITLIST *cykhitlist, int *msamap, int firstpos, STATSMETHOD statsmethod)
{
  int h;
  int ih, jh;

  if (!fp)         return eslOK;
  if (!cykhitlist) return eslOK;

  for (h = 0; h < nhit; h++) cykhitlist->srthit[h] = cykhitlist->hit + h;
  if (nhit > 1) qsort(cykhitlist->srthit, nhit, sizeof(HIT *), (statsmethod == NAIVE)? hit_sorted_by_score:hit_sorted_by_eval);

  fprintf(fp, "# in_cyk in_given   left_pos       right_pos      score           E-value    substitutions      power\n");
  fprintf(fp, "#----------------------------------------------------------------------------------------------------------------\n");
  if (nhit == 0) fprintf(fp, "no significant pairs\n");

  for (h = 0; h < nhit; h ++) {
    ih = hitlist->srthit[h]->i;
    jh = hitlist->srthit[h]->j;
    
    if (hitlist->srthit[h]->bptype == WWc) {
      if (cykhitlist->srthit[h]->bptype == WWc) 
	fprintf(fp, "*\t*\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else if (cykhitlist->srthit[h]->is_compatible)
	fprintf(fp, "~\t*\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else 
	fprintf(fp, " \t*\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
   }
    
    else if (hitlist->srthit[h]->bptype < STACKED) { 
      if (cykhitlist->srthit[h]->bptype == WWc) 
	fprintf(fp, "*\t**\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else if (cykhitlist->srthit[h]->is_compatible)
	fprintf(fp, "~\t**\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else 
	fprintf(fp, " \t**\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
  }
    
    else if (hitlist->srthit[h]->bptype < BPNONE && hitlist->srthit[h]->is_compatible) {
           if (cykhitlist->srthit[h]->bptype == WWc) 
	fprintf(fp, "*\tc~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else if (cykhitlist->srthit[h]->is_compatible)
	fprintf(fp, "~\tc~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else 
	fprintf(fp, " \tc~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
   }
    
    else if (hitlist->srthit[h]->bptype < BPNONE) {
           if (cykhitlist->srthit[h]->bptype == WWc) 
	fprintf(fp, "*\tc\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else if (cykhitlist->srthit[h]->is_compatible)
	fprintf(fp, "~\tc\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else 
	fprintf(fp, " \tc\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
    }
    
    else if (hitlist->srthit[h]->is_compatible) {
           if (cykhitlist->srthit[h]->bptype == WWc) 
	fprintf(fp, "*\t~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else if (cykhitlist->srthit[h]->is_compatible)
	fprintf(fp, "~\t~\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else 
	fprintf(fp, " \t~%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
   }
    
    else {
           if (cykhitlist->srthit[h]->bptype == WWc) 
	fprintf(fp, "*\t \t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else if (cykhitlist->srthit[h]->is_compatible)
	fprintf(fp, "~\t?\t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
      else 
	fprintf(fp, " \t \t%10d\t%10d\t%.5f\t%g\t%lld\t\t%.2f\n", 
		msamap[ih]+firstpos, msamap[jh]+firstpos, hitlist->srthit[h]->sc, hitlist->srthit[h]->Eval, hitlist->srthit[h]->nsubs, hitlist->srthit[h]->power);
    }  
  }

  return eslOK;
}


void
cov_FreeRankList(RANKLIST *ranklist)
{
  if (ranklist == NULL) return;

  if (ranklist->ha)      esl_histogram_Destroy(ranklist->ha);
  if (ranklist->ht)      esl_histogram_Destroy(ranklist->ht);
  if (ranklist->hb)      esl_histogram_Destroy(ranklist->hb);
  if (ranklist->survfit) free(ranklist->survfit);
  free(ranklist);
}

void
cov_FreeHitList(HITLIST *hitlist)
{
  if (hitlist == NULL) return;

  if (hitlist->srthit) free(hitlist->srthit);
  if (hitlist->hit)    free(hitlist->hit);
  free(hitlist);
}

int
cov_SignificantPairs_ZScore(struct mutual_s *mi, int *msamap, int firstpos, int *ct, int verbose, char *errbuf)
{
  double       avgi, avgj;
  double       stdi, stdj;
  double       zscorei, zscorej;
  double       zscore;
  int          i, j, k;
  int          ipair;

  for (i = 0; i < mi->alen; i ++) {
    if (ct[i+1] > 0 && ct[i+1] > i+1) {
      ipair = ct[i+1]-1;

      avgi = 0.0;
      stdi = 0.0;
      
      for (j = 0; j < mi->alen; j ++) {
	if (j != ipair) {   
	  if (j != i) {
	    avgi += mi->COV->mx[i][j];
	    stdi += mi->COV->mx[i][j] * mi->COV->mx[i][j];
	  }
	}
	else {
	  avgj = 0.0;
	  stdj = 0.0;
	  for (k = 0; k < mi->alen; k ++) {
	    if (k != j && k != i) {       
	      avgj += mi->COV->mx[j][k];
	      stdj += mi->COV->mx[j][k] * mi->COV->mx[j][k];
	    }
	  }
	  avgj /= (mi->alen-2);
	  stdj /= (mi->alen-2);
	  stdj -= avgj*avgj;
	  stdj = sqrt(stdj);
	}
      }
      avgi /= (mi->alen-2);
      stdi /= (mi->alen-2);
      stdi -= avgi*avgi;
      stdi = sqrt(stdi);

      zscorei = (mi->COV->mx[i][ipair] - avgi) / stdi;
      zscorej = (mi->COV->mx[i][ipair] - avgj) / stdj;
      zscore  = ESL_MIN(zscorej, zscorej);
      printf("[%d][%d] %f | %f | %f %f | %f %f\n", msamap[i]+firstpos, msamap[ipair]+firstpos, zscore, mi->COV->mx[i][ipair], avgi, stdi, avgj, stdj);
    }
  }  
  return eslOK;
}

int
cov_FisherExactTest(double *ret_pval, int cBP, int cNBP, int BP, int alen)
{
  double pval = 0.0;
  double add, add_fixed;
  double factorial_total, factorial_BP, factorial_NBP;
  double factorial_cnbp,  factorial_ncnbp;
  double factorial_cbp,   factorial_ncbp;
  double factorial_cov,   factorial_non;
  double tol = 1e-3;
  int    NBP;
  int    total;
  int    cov;
  int    non;
  int    cbp, ncbp;
  int    cnbp, ncnbp;
 
  NBP   = alen * (alen - 1) / 2;
  total = BP    + NBP;
  cov   = cBP   + cNBP;
  non   = total - cov;

  esl_stats_LogGamma(total+1, &factorial_total);
  esl_stats_LogGamma(BP+1,    &factorial_BP);
  esl_stats_LogGamma(NBP+1,   &factorial_NBP);
  esl_stats_LogGamma(cov+1,   &factorial_cov);
  esl_stats_LogGamma(non+1,   &factorial_non);

  add_fixed  = factorial_BP + factorial_NBP + factorial_cov + factorial_non - factorial_total;

  for (cnbp = cNBP; cnbp <= NBP; cnbp ++) {
    cbp   = cov - cnbp; if (cbp < 0) break;
    ncbp  = BP  - cbp;
    ncnbp = non - ncbp;
    
    esl_stats_LogGamma(cbp+1,   &factorial_cbp);
    esl_stats_LogGamma(ncbp+1,  &factorial_ncbp);
    esl_stats_LogGamma(cnbp+1,  &factorial_cnbp);
    esl_stats_LogGamma(ncnbp+1, &factorial_ncnbp);

    add = add_fixed - factorial_cbp - factorial_ncbp - factorial_cnbp - factorial_ncnbp;

    pval += exp(add);
  }
  if (pval > 1.0 && pval < 1.0 + tol) pval = 1.0;

  *ret_pval = pval;

  return eslOK;
}


int
cov_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h, double *survfit)
{
  int i;
  uint64_t c = 0;
  double   esum;
  double ai;
  
  /* The observed binned counts:
   */
  if (h->imax >= 0 && h->obs[h->imax] > 1) {
    if (fprintf(fp, "%f\t%g\n", h->xmax, 1.0 / (double) h->Nc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
  }
  for (i = h->imax; i >= h->imin; i--)
    {
      if (h->obs[i] > 0) {
	c   += h->obs[i];
	ai = esl_histogram_Bin2LBound(h, i);
	if (fprintf(fp, "%f\t%g\n", ai, (double) c / (double) h->Nc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
      }
    }
  if (fprintf(fp, "&\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");

  /* The survival fit
   */
  if (survfit != NULL) 
    {
      for (i = 2*h->nb-1; i >= 0; i--)
	{
	  if (survfit[i] > 0.) {
	    esum = survfit[i];
	    ai = esl_histogram_Bin2LBound(h, i);
	    if (fprintf(fp, "%f\t%g\n", ai, esum) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
	  }
	}
      if (fprintf(fp, "&\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
    }
  return eslOK;
}

int
cov_histogram_SetSurvFitTail(ESL_HISTOGRAM *h, double **ret_survfit, double pmass, double (*surv)(double x, void *params), void *params)
{
  double *survfit = NULL;
  int     status;
  int     b;
  double  ai, bi;

  ESL_ALLOC(survfit, sizeof(double) * 2*h->nb);
  esl_vec_DSet(survfit, 2*h->nb, 0.);

  for (b = h->cmin; b < 2*h->nb; b++)
    {
      ai = esl_histogram_Bin2LBound(h, b);
      bi = esl_histogram_Bin2UBound(h, b);
      survfit[b] = pmass * (*surv)(bi, params);
    }

  *ret_survfit = survfit;
  return eslOK;

 ERROR:
  return status;
}

int 
cov_WriteHistogram(struct data_s *data, char *gnuplot, char *covhisfile, char *covqqfile, SAMPLESIZE samplesize, RANKLIST *ranklist, char *title)
{
  FILE     *fp = NULL;
  RANKLIST *ranklist_null = data->ranklist_null;
  char     *errbuf = data->errbuf;
  int       ignorebps;
  int       status;

  if (ranklist == NULL) return eslOK;

  if (covhisfile) {
    if ((fp = fopen(covhisfile, "w")) == NULL) ESL_XFAIL(eslFAIL, errbuf, "could not open covhisfile %s\n", covhisfile);
    /* write the survival for the given alignment (base pairs) */
    cov_histogram_PlotSurvival(fp, ranklist->hb, NULL);
    /* write the survival for the given alignment (not pairs) */
    cov_histogram_PlotSurvival(fp, ranklist->ht, NULL);
    /* write the survival for the null alignments */
    if (ranklist_null) cov_histogram_PlotSurvival(fp, ranklist_null->ha, ranklist_null->survfit);
    fclose(fp);

    ignorebps = (samplesize == SAMPLE_ALL)? TRUE : data->ignorebps;
    if (!data->nofigures) {
      status = cov_PlotHistogramSurvival(data, gnuplot, covhisfile, ranklist, title, FALSE, ignorebps);
      if (status != eslOK) goto ERROR;
      status = cov_PlotHistogramSurvival(data, gnuplot, covhisfile, ranklist, title, TRUE,  ignorebps);
      if (status != eslOK) goto ERROR;
    }
  }

#if 0
  // the qq plots. Not as useful as the histogram plots above. 
  // commented to avoid extra time
  status = cov_PlotHistogramQQ(data, gnuplot, covqqfile, ranklist, title, FALSE);
  if (status != eslOK) goto ERROR;
  status = cov_PlotHistogramQQ(data, gnuplot, covqqfile, ranklist, title, TRUE);
  if (status != eslOK) goto ERROR;
#endif

  return eslOK;

 ERROR:
  return status;

}

int 
cov_NullFitExponential(ESL_HISTOGRAM *h, double **ret_survfit, double pmass, double *ret_newmass, double *ret_mu, double *ret_lambda, int verbose, char *errbuf)
{
  double         ep[2];  	/* estimated mu, lambda  */
  double         newmass;
  int            status;

  /* set the tail by mass */
  status = esl_histogram_SetTailByMass(h, pmass, &newmass);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not set TailByMass");

  /* exponential fit to tail */
  status = esl_exp_FitCompleteBinned(h, &ep[0], &ep[1]);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not do exponential fit");
  
  /* the expected survival */
  if (!isinf(ep[1])) {
    cov_histogram_SetSurvFitTail(h, ret_survfit, newmass, &esl_exp_generic_surv, ep);
  }
  
  *ret_mu      = ep[0];
  *ret_lambda  = ep[1];
  *ret_newmass = newmass;
  return eslOK;

 ERROR:
  return status;
}

int 
cov_NullFitGamma(ESL_HISTOGRAM *h, double **ret_survfit, double pmass, double *ret_newmass, double *ret_mu, double *ret_lambda, double *ret_k, int verbose, char *errbuf)
{
  double         ep[3];  	/* ep[0] = mu; ep[1]=lambda=1/2; ep[2] = tau = k/2 */
  double         newmass;
  int            status;

  /* set the tail by mass */
  status = esl_histogram_SetTailByMass(h, pmass, &newmass);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not set TailByMass");

  /* Gamma fit to tail */
  status = esl_gam_FitCompleteBinned(h, &ep[0], &ep[1], &ep[2]);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "could not do a Gamma fit");
  if (verbose) printf("GammaFit mu = %f lambda = %f tau = %f (mass %f)\n", ep[0], ep[1], ep[2], newmass);

  /* the expected survival, if we do have a fit */
  if (!isinf(ep[1])) {
    cov_histogram_SetSurvFitTail(h, ret_survfit, newmass, &esl_gam_generic_surv, ep);
  }
 
  *ret_mu     = ep[0];
  *ret_lambda = ep[1];
  *ret_k      = ep[2];
  *ret_newmass = newmass;
  return eslOK;

 ERROR:
  return status;
}



int 
cov_PlotHistogramSurvival(struct data_s *data, char *gnuplot, char *covhisfile, RANKLIST *ranklist, char *title, int dosvg, int ignorebps)
{
  FILE     *pipe;
  RANKLIST *ranklist_null = data->ranklist_null;
  RANKLIST *ranklist_aux  = data->ranklist_aux;
  char     *filename = NULL;
  char     *outplot = NULL;
  char     *key0 = NULL;
  char     *key1 = NULL;
  char     *key2 = NULL;
  char     *key3 = NULL;
  char     *key4 = NULL;
  char     *key  = NULL;
  double    minphi;
  double    minmass = 0.005;
  int       pointype;
  double    pointintbox;
  int       linew;
  double    pointsize;
  double    xmin, xmax;
  double    ymin, ymax;
  double    y2min, y2max;
  double    posx, posy;
  double    incx, incy;
  double    expsurv = data->thresh->val;
  double    extra;
  int       has_bpairs = FALSE;
  int       linespoints;
  int       status;
  
  if (gnuplot    == NULL) return eslOK;
  if (covhisfile == NULL) return eslOK;
  if (ranklist   == NULL) return eslOK;
  
  if (data->mode == GIVSS && data->nbpairs     > 0) has_bpairs = TRUE;
  if (data->mode == CYKSS && data->nbpairs_cyk > 0) has_bpairs = TRUE;

  esl_FileTail(covhisfile, FALSE, &filename);
  
  esl_sprintf(&key0, "all pairs");
  esl_sprintf(&key1, "not proposed pairs");
  esl_sprintf(&key2, "proposed pairs");
  esl_sprintf(&key3, "null distribution");
  if (ranklist_aux) esl_sprintf(&key4, "null-null distribution");
  
  /* for plotting */
  esl_histogram_SetTailByMass(ranklist->ha, minmass, NULL);
  minphi = ranklist->ha->phi; 
  
  /* do the plotting */
  pipe = popen(gnuplot, "w");
  if (dosvg) {
    esl_sprintf(&outplot, "%s.svg", covhisfile);
    fprintf(pipe, "set terminal svg font 'Arial,12'\n");
    pointype    = 66;
    pointsize   = 0.6;
    pointintbox = 0.4;
    linew       = 1;
  }
  else {
    esl_sprintf(&outplot, "%s.ps", covhisfile);
    fprintf(pipe, "set terminal postscript color 14\n");
    pointype    = 65;
    pointsize   = 0.7;
    pointintbox = 1.0;
    linew       = 2;
  }
  
  fprintf(pipe, "set output '%s'\n", outplot);
  fprintf(pipe, "set title '%s' noenhanced\n", title);
  
  fprintf(pipe, "set style line 1   lt 1 lc rgb 'red'        pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 2   lt 1 lc rgb '#DC143C'    pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 3   lt 1 lc rgb 'blue'       pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 4   lt 1 lc rgb '#00008B'    pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 5   lt 1 lc rgb 'black'      pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 11  lt 1 lc rgb 'red'        pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 22  lt 1 lc rgb '#DC143C'    pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 33  lt 1 lc rgb 'blue'       pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 44  lt 1 lc rgb '#00008B'    pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 55  lt 1 lc rgb 'black'      pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 111 dashtype 3 lc rgb 'red'  pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 333 dashtype 3 lc rgb 'blue' pt 7 lw %d ps %f\n", linew, pointsize);
    
  // plot evalue
  fprintf(pipe, "set multiplot layout 1,1\n");  
  fprintf(pipe, "set lmargin 5\n");  
  if (dosvg) fprintf(pipe, "set lmargin 10\n");
  fprintf(pipe, "set bmargin 1\n");  
  if (dosvg) fprintf(pipe, "set bmargin 4\n");
  fprintf(pipe, "set rmargin 6\n");
  if (dosvg) fprintf(pipe, "set rmargin 9\n");
  fprintf(pipe, "set tmargin 1\n");  
  if (dosvg) fprintf(pipe, "set tmargin 4\n");
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");  
  fprintf(pipe, "set logscale y\n");
  fprintf(pipe, "set logscale y2\n");
  fprintf(pipe, "set format y   '%s'\n", "10^{%T}");
  fprintf(pipe, "set format y2  '%s' \n", "10^{%T}");
  fprintf(pipe, "set ytics nomirror \n");
  fprintf(pipe, "set ylabel  'Distribution of pairs with covariation score > t'\n");
  fprintf(pipe, "set xlabel 'covariation score t'\n");
  if (!ignorebps && has_bpairs) {
    fprintf(pipe, "set y2tics textcolor ls 3 \n");
    fprintf(pipe, "set y2label 'Expected # proposed pairs with score > t' offset -1,0 textcolor ls 3 \n");
  }
  
  // the ymax and xmax values
  xmax = ESL_MAX(ranklist->hb->xmax, ranklist->ht->xmax);
  if (ranklist_null) {
    xmax = ESL_MAX(xmax, evalue2cov(expsurv, ranklist->ht->Nc, ranklist_null->ha, ranklist_null->survfit));
    xmax = ESL_MAX(xmax, evalue2cov(expsurv, ranklist->hb->Nc, ranklist_null->ha, ranklist_null->survfit));
  }
  xmin = (ranklist_null)?
    ESL_MIN(ranklist->ht->xmin, evalue2cov(expsurv,ranklist->ht->Nc, ranklist_null->ha, ranklist_null->survfit)) : ranklist->ht->xmin;
  
  extra = (xmax-xmin)/10.;
  xmax += extra;
  xmin -= extra/2.;

  ymax = (ranklist_null)? cov2evalue(xmin, 1, ranklist_null->ha, ranklist_null->survfit) + 0.7 : 1.7;
  ymin = (ranklist_null)? cov2evalue(xmax, 1, ranklist_null->ha, ranklist_null->survfit)       : 1.0/(double)ranklist->ha->Nc;
  ymin = ESL_MIN(ymin, expsurv/(double)ranklist->ht->Nc);
  ymin = ESL_MIN(ymin, expsurv/(double)ranklist->hb->Nc);
  ymin *= exp(-0.1*(fabs(log(ymax) - log(ymin))));

  if (ymin <= 1e-15) ymin = 1e-15;
  if (ymin == ymax) ymin = 1e-5;

  incx = (xmax-xmin)/10.;
  incy = fabs(log(ymax)-log(ymin))/15.;
  xmax += incx;
  
  y2min = ymin * ranklist->ht->Nc;
  y2max = ymax * ranklist->ht->Nc;
  fprintf(pipe, "set xrange   [%f:%f]\n", xmin, xmax);
  fprintf(pipe, "set yrange   [%g:%g]\n", ymin, ymax);
  if (!ignorebps && has_bpairs) {
    fprintf(pipe, "set y2range  [%g:%g]\n", (double)ranklist->hb->Nc*ymin, (double)ranklist->hb->Nc*ymax);
  }
  
  posx = xmin + 2.2*incx;

  if (ranklist_null) { 
    esl_sprintf(&key, "E %.3f", expsurv);

    if (ignorebps) {
      cov_plot_lineatexpcov  (pipe, data, expsurv, ranklist->ha->Nc, ranklist_null->ha, ranklist_null->survfit, ranklist->ha, "x1y1",
			      key, ymax, ymin, xmax, xmin, 3, 333);
     }
    else {
      cov_plot_lineatexpcov  (pipe, data, expsurv, ranklist->ht->Nc, ranklist_null->ha, ranklist_null->survfit, ranklist->ht, "x1y1",
			      key, ymax, ymin, xmax, xmin, 1, 111);
      if (has_bpairs)
	cov_plot_lineatexpcov(pipe, data, expsurv, ranklist->hb->Nc, ranklist_null->ha, ranklist_null->survfit, ranklist->hb, "x1y1",
			      key, ymax, ymin, xmax, xmin, 3, 333);
    }
    
    linespoints = FALSE;
    posy = ymin*exp(1.0*incy);

    status = cov_histogram_plotexpectsurv(pipe, ranklist_null->ha->Nc, ranklist_null->ha, ranklist_null->survfit, key3, "x1y1",
					  FALSE, posx, posy, FALSE, linespoints, 55, 5);
    if (status != eslOK) goto ERROR;
    linespoints = TRUE;
  }

  if (ignorebps) {
    posy = ymin*exp(2.0*incy);
    status = cov_histogram_plotexpectsurv(pipe, ranklist->ha->Nc, ranklist->ha, NULL, key0, "x1y1", FALSE, posx, posy,
					  FALSE, linespoints, 33, 3);
    if (status != eslOK) goto ERROR;
    cov_plot_extra_yaxis(pipe, (double)ranklist->ht->Nc*ymax, (double)ranklist->ht->Nc*ymin, -4.0, "Expected # pairs with score > t", 3);
  }
  else {
    posy = ymin*exp(2.0*incy);
    status = cov_histogram_plotexpectsurv  (pipe, ranklist->ht->Nc, ranklist->ht, NULL, key1, "x1y1", FALSE, posx, posy,
					    FALSE, linespoints, 11, 1);
    if (status != eslOK) goto ERROR;
    
    if (has_bpairs) { // the distribution of base-pairs pairs
      posy = ymin*exp(3.0*incy);
      status = cov_histogram_plotexpectsurv(pipe, ranklist->hb->Nc, ranklist->hb, NULL, key2, "x1y1", FALSE, posx, posy,
					    FALSE, linespoints, 33, 3);
      if (status != eslOK) goto ERROR;
    }
    cov_plot_extra_yaxis(pipe, (double)ranklist->ht->Nc*ymax, (double)ranklist->ht->Nc*ymin, -12.0, "Expected # not proposed pairs with score > t", 1);
  }
  
  pclose(pipe);
  
  free(key0);
  free(key1);
  free(key2);
  free(key3);
  free(key4);
  free(key);
  free(outplot);
  free(filename);
  return eslOK;
  
 ERROR:
  if (key0) free(key0);
  if (key1) free(key1);
  if (key2) free(key2);
  if (key3) free(key3);
  if (key4) free(key4);
  if (key)  free(key);
  if (outplot) free(outplot);
  if (filename) free(filename);
  pclose(pipe);
  return status;
}

int 
cov_PlotHistogramQQ(struct data_s *data, char *gnuplot, char *covhisfile, RANKLIST *ranklist, char *title, int dosvg)
{
  FILE     *pipe;
  RANKLIST *ranklist_null = data->ranklist_null;
  char     *filename = NULL;
  char     *outplot = NULL;
  char     *key1 = NULL;
  char     *key2 = NULL;
  double    minphi;
  double    minmass = 0.005;
  int       pointype;
  double    pointintbox;
  int       linew;
  double    pointsize;
  double    max;
  double    en, eo;
  double    cov;
  int       linespoints = TRUE;
  int       nsample = 0;
  int       status;

  if (gnuplot    == NULL) return eslOK;
  if (covhisfile == NULL) return eslOK;
  if (ranklist   == NULL) return eslOK;

  esl_FileTail(covhisfile, FALSE, &filename);

  esl_sprintf(&key1, "all pairs");
  esl_sprintf(&key2, "not base pairs");

  /* for plotting */
  esl_histogram_SetTailByMass(ranklist->ha, minmass, NULL);
  minphi = ranklist->ha->phi; 
 
  /* do the plotting */
  pipe = popen(gnuplot, "w");
  if (dosvg) {
    esl_sprintf(&outplot, "%s.svg", covhisfile);
    fprintf(pipe, "set terminal svg font 'Arial,12'\n");
    pointype    = 71;
    pointsize   = 0.6;
    pointintbox = 0.4;
    linew       = 1;
  }
  else {
    esl_sprintf(&outplot, "%s.ps", covhisfile);
    fprintf(pipe, "set terminal postscript color 14\n");
    pointype    = 65;
    pointsize   = 0.7;
    pointintbox = 1.0;
    linew       = 2;
  }

  fprintf(pipe, "set output '%s'\n", outplot);
  fprintf(pipe, "set title '%s' noenhanced \n", title);

  fprintf(pipe, "set style line 1   lt 1 lc rgb 'grey'      pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 2   lt 1 lc rgb 'brown'     pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 3   lt 1 lc rgb 'cyan'      pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 4   lt 1 lc rgb 'red'       pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 5   lt 1 lc rgb 'orange'    pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 6   lt 1 lc rgb 'turquoise' pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 7   lt 1 lc rgb 'black'     pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 8   lt 1 lc rgb 'green'     pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 9   lt 1 lc rgb 'blue'      pt 7 lw %d ps %f\n", linew, pointsize);
  fprintf(pipe, "set style line 11  lt 1 lc rgb 'grey'      pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 22  lt 1 lc rgb 'brown'     pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 33  lt 1 lc rgb 'cyan'      pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 44  lt 1 lc rgb 'red'       pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 55  lt 1 lc rgb 'orange'    pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 66  lt 1 lc rgb 'turquoise' pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 77  lt 1 lc rgb 'black'     pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 88  lt 1 lc rgb 'green'     pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);
  fprintf(pipe, "set style line 99  lt 1 lc rgb 'blue'      pt %d pi -1  lw %d ps %f \nset pointintervalbox %f\n", pointype, linew, pointsize, pointintbox);

  // plot qq
  fprintf(pipe, "set multiplot\n");  
  fprintf(pipe, "set ylabel 'Observed #pairs'\n");
  fprintf(pipe, "set xlabel 'Expected #pairs (E-value)'\n");

  // the max values
  max = 50.;
  cov = (ranklist_null)? ESL_MAX(ranklist->ha->xmax,ranklist_null->ha->xmax) : ranklist->ha->xmax;
  eo = cov2evalue(cov, nsample, ranklist->ha, NULL);
  en = cov2evalue(cov, nsample, ranklist_null->ha, ranklist_null->survfit);

  fprintf(pipe, "set logscale x\n");
  //fprintf(pipe, "set logscale y\n");
  fprintf(pipe, "set yrange [%g:%f]\n", 0.0, max);
  fprintf(pipe, "set xrange [%g:%f]\n", en, 2.0*max);

  fprintf(pipe, "plot x title 'exp=obs' ls 77\n");
  if (ranklist_null) {
    status = cov_histogram_plotqq(pipe, data, ranklist->ha, ranklist_null->ha, key1, FALSE, linespoints, 99, 2);
    if (status != eslOK) goto ERROR;
#if 0
    status = cov_histogram_plotqq(pipe, data, ranklist->ht, ranklist_null->ha, key2, FALSE, linespoints, 44, 2);
    if (status != eslOK) goto ERROR;
#endif
  }
  pclose(pipe);
  
  free(key1);
  free(key2);
  free(outplot);
  free(filename);
  return eslOK;

 ERROR:
  if (key1) free(key1);
  if (key2) free(key2);
  if (outplot) free(outplot);
  if (filename) free(filename);
  if (pipe) pclose(pipe);
  return status;
}






int
cov_ranklist_Bin2Bin(int b, ESL_HISTOGRAM *h, ESL_HISTOGRAM *new, int *ret_newb)
{
  double x;
  int    newb;
  int    status;
  
  x = esl_histogram_Bin2LBound(h, b);
  if (! isfinite(x)) ESL_XEXCEPTION(eslERANGE, "value added to histogram is not finite");

  x = round( ((x - new->bmin) / new->w)); 

  /* x is now the bin number as a double, which we will convert to
   * int. Because x is a double (64-bit), we know all ints are exactly
   * represented.  Check for under/overflow before conversion.
   */
  if (x < (double) INT_MIN || x > (double) INT_MAX) 
    ESL_XEXCEPTION(eslERANGE, "value %f isn't going to fit in histogram", x);
  
  newb = (int) x;
  if (newb > new->nb) 
    ESL_XEXCEPTION(eslERANGE, "bin value %d isn't going to fit in histogram bin %d ", newb, new->nb);
  
  *ret_newb = newb;
  return eslOK;

 ERROR:
  *ret_newb = -1;
  return status;
}


/*---------------- internal functions --------------------- */



static double
cov2evalue(double cov, int Nc, ESL_HISTOGRAM *h, double *survfit)
{
  double eval = +eslINFINITY;
  int    icov;
  int    c = 0;
  int    i;

  esl_histogram_Score2Bin(h, cov, &icov);

  /* use the fit if possible */
  if      (survfit && icov >= 2*h->nb-1)  eval = survfit[2*h->nb-1] * (double)Nc;
  else if (survfit &&  cov >= h->phi)     eval = survfit[icov+1]    * (double)Nc;
  else { /* otherwise, the sampled distribution  */
    if (cov >= h->xmax) { return (double)Nc / (double)h->Nc; }

    if (icov <= h->imax) {      
      if (icov <  h->imin)   icov = h->imin;
      if (icov >= h->imax-1) eval = (double)Nc / (double)h->Nc;
      else {
	for (i = h->imax; i >= icov; i--) c += h->obs[i];
	eval = (double)c * (double)Nc / (double)h->Nc;
      }
    }
    else {
      printf("cannot find evalue for covariation %f\n", cov); exit(1);
    }
  }

  return eval;
}

static double
evalue2cov(double eval, int Nc, ESL_HISTOGRAM *h, double *survfit)
{
  double cov = -eslINFINITY;
  double exp = 0.0;
  int    min_nobs = 10;
  int    c = 0;
  int    i;
  int    b = h->cmin-1; 
  
  /* use the fit if possible */
  if (h->No >= min_nobs && survfit) {
    for (b = 2*h->nb-1; b >= h->cmin; b--) {
      exp = survfit[b];
      if (exp * (double)Nc >= eval) break;
    }
    cov = esl_histogram_Bin2LBound(h, b+1);
  }

  /* otherwise, use the sampled distribution */
  if (b == h->cmin-1) {
    for (i = h->imax; i >= h->imin; i--) {
      c += h->obs[i];
      if ((double)c * (double)Nc / (double)h->Nc > eval) break;
    }    
    cov = esl_histogram_Bin2LBound(h, i+1);
  }
  
  return cov;
}

static double
cov_histogram_pmass(ESL_HISTOGRAM *h, double target_pmass, double target_fracfit)
{
  int       i;
  int       tp = 0;   // total number of points in histogram with observations
  int       nfit = 0;
  uint64_t  c = 0;
  double    pmass;

  for (i = h->imax; i >= h->imin; i--) { if (h->obs[i] > 0) tp ++; }
  i = h->imax;
	 
  for (i = h->imax; i >= h->imin; i--)
    {
      c += h->obs[i];
      if (h->obs[i] > 0) {
	nfit ++;
	pmass = (double)c / (double)h->Nc;
	if ((double)nfit/(double)tp >= target_fracfit || pmass >= target_pmass) break;
      }
    }
  return pmass;
}


static int
cov_histogram_plotdensity(FILE *pipe, ESL_HISTOGRAM *h, double *survfit, char *key, double posx, double posy, int logval, int style1, int style2)
{
  int       i;
  uint64_t  obs;
  double    exp;
  double    ai;
 
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "set label 1 at %f,%f '%s' center tc ls %d\n", posx, posy, key, style1);
  fprintf(pipe, "plot '-' using 1:2 with linespoints ls %d \n", style1);

  if (h->obs[h->imax] > 1) 
    if (fprintf(pipe, "%f\t%f\n", 
		h->xmax, (logval)? -log((double)h->Nc) : 1.0/(double) h->Nc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
  for (i = h->imax; i >= h->imin; i--)
    {
      obs = h->obs[i];

      if (obs > 0) {
	ai = esl_histogram_Bin2LBound(h, i);
	if (fprintf(pipe, "%f\t%g\n", 
		    ai, (logval)? log((double)obs)-log((double)h->Nc) : (double)obs/(double) h->Nc) < 0) 
	  ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
      }
    }
  fprintf(pipe, "e\n");
  
  /* The expected binned counts:
   */
  if (survfit != NULL) 
    {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");
      fprintf(pipe, "set key off\n");
      fprintf(pipe, "plot '-' using 1:2 with lines ls %d \n", style2);
      
      for (i = h->nb-1; i > 0; i--)
	{
	  exp = survfit[i]-survfit[i-1];        /* some worry about 1+eps=1 problem here */

	  if (exp > 0.) { 
	    ai = esl_histogram_Bin2LBound(h, i);
	    if (fprintf(pipe, "%f\t%g\n", 
			ai, (logval)? log(exp)-log((double)h->Nc) : exp/(double) h->Nc) < 0) 
	      ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
	  }
	}
      fprintf(pipe, "e\n"); 
    }
  
  return eslOK;
}
static int
cov_histogram_plotsurvival(FILE *pipe, ESL_HISTOGRAM *h, double *survfit, char *key, double posx, double posy, int logval, int style1, int style2)
{
  int       i;
  uint64_t  c = 0;
  double    esum;
  double    ai;
 
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "set label 1 at %f,%f '%s' center tc ls %d\n", posx, posy, key, style1);
  fprintf(pipe, "plot '-' using 1:2 with linespoints ls %d \n", style1);

  for (i = h->imax; i >= h->imin; i--)
    {
      c += h->obs[i];

      if (h->obs[i] > 0) {
	ai = esl_histogram_Bin2LBound(h, i);
	if (fprintf(pipe, "%f\t%g\n", 
		    ai, (logval)? log((double)c)-log((double)h->Nc) : (double)c/(double) h->Nc) < 0) 
	  ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
      }
    }
  fprintf(pipe, "e\n");
  
  /* The survival fit:
   */
  if (survfit != NULL) 
    {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");
      fprintf(pipe, "set key off\n");
      fprintf(pipe, "plot '-' using 1:2 with lines ls %d \n", style2);
      
      esum = 0.;
      for (i = h->nb-1; i >= 0; i--)
	{
	  esum = survfit[i];    
	  
	  if (esum > 0.) { 
	    ai = esl_histogram_Bin2LBound(h, i);
	    if (fprintf(pipe, "%f\t%g\n", 
			ai, (logval)? log(esum) : esum) < 0) 
	      ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
	  }
	}
      fprintf(pipe, "e\n"); 
    }
  
  return eslOK;
}

static int
cov_histogram_plotqq(FILE *pipe, struct data_s *data, ESL_HISTOGRAM *h1, ESL_HISTOGRAM *h2, char *key, int logval, int linespoints, int style1, int style2)
{
  int       i;
  uint64_t  c = 0;
  double    eval = 0.;
  double    ai;
 
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "plot '-' using 1:2 title '%s' with linespoints ls %d \n", key, style1);

  for (i = h1->imax; i >= h1->imin; i--)
    {
      c += h1->obs[i];
  
      if (h1->obs[i] > 0) {
	ai = esl_histogram_Bin2LBound(h1, i);
	eval = cov2evalue(ai, h1->Nc, h2, NULL);
	if (fprintf(pipe, "%g\t%g\n", 
		    (logval)? log(eval)      : eval,
		    (logval)? log((double)c) : (double)c) < 0) 
	  ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
     }
    }
  fprintf(pipe, "e\n");
   
  return eslOK;
}

static int
cov_histogram_plotexpectsurv(FILE *pipe, int Nc, ESL_HISTOGRAM *h, double *survfit, char *key, char *axes, int addtics, double posx, double posy, int logval, 
			     int linespoints, int style1, int style2)
{
  int       i;
  uint64_t  c = 0;
  double    esum;
  double    ai;

  if (posy <= 0.) posy = 1e-5;
  
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "set label 1 at %f,%g '%s' center tc ls %d\n", posx, posy, key, style1);
  if (addtics) {
    if (linespoints) fprintf(pipe, "plot '-' using 1:2:yticlabels(3) axes %s with linespoints ls %d \n", axes, style1);
    else             fprintf(pipe, "plot '-' using 1:2:yticlabels(3) axes %s with points ls %d \n", axes, style1);
  }
  else {
    if (linespoints) fprintf(pipe, "plot '-' using 1:2 axes %s with linespoints ls %d \n", axes, style1);
    else             fprintf(pipe, "plot '-' using 1:2 axes %s with points ls %d \n", axes, style1);
  }

  for (i = h->imax; i >= h->imin; i--)
    {
      c += h->obs[i];
      if (h->obs[i] > 0) {
	ai = esl_histogram_Bin2LBound(h, i);
	if (fprintf(pipe, "%g\t%g\t%g\n", 
		    ai,
		    (logval)? log((double)c)                   - log((double)h->Nc) : (double)c              / (double)h->Nc,
		    (logval)? log((double)c) + log((double)Nc) - log((double)h->Nc) : (double)c * (double)Nc / (double)h->Nc) < 0) 
	  ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
      }
    }
  fprintf(pipe, "e\n");
  
  /* The survplot:
   */
  if (survfit != NULL) 
    {
      fprintf(pipe, "set size 1,1\n");
      fprintf(pipe, "set origin 0,0\n");
      fprintf(pipe, "set key off\n");
      fprintf(pipe, "plot '-' using 1:2 with lines ls %d \n", style2);
      
      for (i = 2*h->nb-1; i >= 0; i--)
	{
	  esum = survfit[i];
	  if (esum > 0.) { 
	    ai = esl_histogram_Bin2LBound(h, i);
	    if (fprintf(pipe, "%g\t%g\t%g\n", 
			ai,
			(logval)? log(esum)                   : esum,
			(logval)? log(esum) + log((double)Nc) : esum * (double)Nc) < 0) 
	      ESL_EXCEPTION_SYS(eslEWRITE, "histogram survival plot write failed");
	  }
	}
      fprintf(pipe, "e\n"); 
    }
  
  return eslOK;
}

static int
cov_plot_lineatexpcov(FILE *pipe, struct data_s *data, double expsurv, int Nc, ESL_HISTOGRAM *h, double *survfit, ESL_HISTOGRAM *h2, char *axes, char *key, 
		      double ymax, double ymin, double xmax, double xmin, int style1, int style2)
{
  double cov;
  double eval;
  double eval2;
  double posx, posy;
  double lenx, leny;
  double ex, ey, dy;
 
  cov   = evalue2cov(expsurv, Nc, h, survfit); if (cov <= -eslINFINITY) return eslOK;
  eval  = cov2evalue(cov, 1, h, survfit);
  eval2 = cov2evalue((cov<h2->xmax)?cov:h2->xmax, 1, h2, NULL);

  lenx = xmax - cov;
  leny = fabs(log(eval) - log(eval2));

  ex = lenx/25.;
  ey = exp(-leny/2.);
  dy = exp(+leny/4.);

  posx = cov + 1.4*ex;
  posy = eval2 * exp(-8.*ey);

  /* The vertical line */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "set label 1 at %g,%g '%s' center tc ls %d\n", posx, posy, key, style1);
  fprintf(pipe, "plot '-' using 1:2 axes %s with lines ls %d \n", axes, style1);
  fprintf(pipe, "%g\t%g\n", cov, (eval*ey >ymin)? eval*ey :ymin);
  fprintf(pipe, "%g\t%g\n", cov, (eval2*dy<ymax)? eval2*dy:ymax);
  fprintf(pipe, "e\n");

  /* The horizontal line */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set key off\n");
  fprintf(pipe, "plot '-' using 1:2 axes %s with lines ls %d \n", axes, style2);
  fprintf(pipe, "%g\t%g\n", (cov-ex > xmin)? cov-ex:xmin, eval);
  fprintf(pipe, "%g\t%g\n", xmax,                         eval);
  fprintf(pipe, "e\n");
  
  return eslOK;
}

static int
cov_plot_extra_yaxis(FILE *pipe, double ymax, double ymin, double xoff, char *label, int style)
{
  /* The observed binned counts:
   */
  fprintf(pipe, "set size 1,1\n");
  fprintf(pipe, "set origin 0,0\n");
  fprintf(pipe, "set noxtic\n");
  fprintf(pipe, "set noxlabel\n");
  fprintf(pipe, "set y2tics offset -5.0,0 nomirror textcolor ls %d \n", style);
  fprintf(pipe, "set y2label '%s' offset %f,0 textcolor ls %d\n", label, xoff, style);
  fprintf(pipe, "set y2range [%g:%g]\n", ymin, ymax);
  fprintf(pipe, "plot -1 notitle \n");
   
  return eslOK;
}

