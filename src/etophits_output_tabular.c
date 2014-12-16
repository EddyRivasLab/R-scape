/* Tabular (parsable) output of pipeline results.
 */
#include "p7_config.h"

#include <stdio.h>
#include <string.h>

#include "hmmer.h"
#include "easel.h"

#include "etophits_output_tabular.h"

/* Function:  p7_tophits_TabularTargets()
 * Synopsis:  Output parsable table of per-sequence hits.
 *
 * Purpose:   Output a parseable table of reportable per-sequence hits
 *            in sorted tophits list <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>, using final pipeline
 *            accounting stored in <pli>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEWRITE> if a write to <ofp> fails; for example, if
 *            the disk fills up.
 */
int
p7_etophits_TabularTargets(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header)
{
  int qnamew = ESL_MAX(20, strlen(qname));
  int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int qaccw  = ((qacc != NULL) ? ESL_MAX(10, strlen(qacc)) : 10);
  int taccw  = ESL_MAX(10, p7_tophits_GetMaxAccessionLength(th));
  int posw   = (pli->long_targets ? ESL_MAX(7, p7_tophits_GetMaxPositionLength(th)) : 0);
  int h,d;

  if (show_header)
  {
      if (pli->long_targets) 
      {
        if (fprintf(ofp, "#%-*s %-*s %-*s %-*s %s %s %*s %*s %*s %*s %*s %6s %9s %6s %5s %s\n",
          tnamew-1, " target name",        taccw, "accession",  qnamew, "query name",           qaccw, "accession", "hmmfrom", "hmm to", posw, "alifrom", posw, "ali to", posw, "envfrom", posw, "env to", posw, ( pli->mode == p7_SCAN_MODELS ? "modlen" : "sq len" ), "strand", "  E-value", " score", " bias", "description of target") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        if (fprintf(ofp, "#%*s %*s %*s %*s %s %s %*s %*s %*s %*s %*s %6s %9s %6s %5s %s\n",
          tnamew-1, "-------------------", taccw, "----------", qnamew, "--------------------", qaccw, "----------", "-------", "-------", posw, "-------", posw, "-------",  posw, "-------", posw, "-------", posw, "-------", "------", "---------", "------", "-----", "---------------------") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-per-sequence hit list: write failed");
      }
      else
      {
        if (fprintf(ofp, "#%*s %33s %22s %33s\n", tnamew+qnamew+taccw+qaccw+2, "", "-------- full sequence --------", "--- best 1 domain ----", "--- domain number estimation ----") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        if (fprintf(ofp, "#%-*s %-*s %-*s %-*s %9s %6s %5s %5s %5s  %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
          tnamew-1, " target name",        taccw, "accession",  qnamew, "query name",           qaccw, "accession",  "  E-value", " score", " bias", " fwd_score", " spfwd_score", " time", "  E-value", " score", " bias", "exp", "reg", "clu", " ov", "env", "dom", "rep", "inc", "description of target") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        if (fprintf(ofp, "#%*s %*s %*s %*s %9s %6s %5s %5s %5s %9s %6s %5s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
          tnamew-1, "-------------------", taccw, "----------", qnamew, "--------------------", qaccw, "----------", "---------", "------", "-----", "-----", "-----", "-----", "---------", "------", "-----", "---", "---", "---", "---", "---", "---", "---", "---", "---------------------") < 0)
          ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
      }
  }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)    
    {
        d    = th->hit[h]->best_domain;
        if (pli->long_targets) 
        {
            if (fprintf(ofp, "%-*s %-*s %-*s %-*s %7d %7d %*d %*d %*d %*d %*" PRId64 " %6s %9.2g %6.1f %5.1f  %s\n",
                tnamew, th->hit[h]->name,
                taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
                qnamew, qname,
                qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
                th->hit[h]->dcl[d].ad->hmmfrom,
                th->hit[h]->dcl[d].ad->hmmto,
                posw, th->hit[h]->dcl[d].ia,
                posw, th->hit[h]->dcl[d].ib,
                posw, th->hit[h]->dcl[d].iae,
                posw, th->hit[h]->dcl[d].ibe,
                posw, th->hit[h]->dcl[0].ad->L,
                (th->hit[h]->dcl[d].ia < th->hit[h]->dcl[d].ib ? "   +  "  :  "   -  "),
                exp(th->hit[h]->lnP),
                th->hit[h]->score,
                th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
                th->hit[h]->desc == NULL ? "-" :  th->hit[h]->desc ) < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        }
        else
        {
                if (fprintf(ofp, "%-*s %-*s %-*s %-*s %9.2g %6.1f %5.1f %9.2g %5.1f %6.1f %5.1f %5.1f %5.1f %5.3f %3d %3d %3d %3d %3d %3d %3d %s\n",
			    tnamew, th->hit[h]->name,
			    taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
			    qnamew, qname,
			    qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
			    exp(th->hit[h]->lnP) * pli->Z,
			    th->hit[h]->score,
			    th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
			    th->hit[h]->fwd_score,
			    th->hit[h]->spfwd_score,
			    th->hit[h]->time,
			    exp(th->hit[h]->dcl[d].lnP) * pli->Z,
			    th->hit[h]->dcl[d].bitscore,
			    th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
			    th->hit[h]->nexpected,
			    0,	/* SRE: FIXME (<nregions> removed now)   */
			    0, 	/* SRE: FIXME (<nclustered> removed now) */
			    th->hit[h]->noverlaps,
			    0,	/* SRE: FIXME (<nenvelopes> removed now) */
			    th->hit[h]->ndom,
			    th->hit[h]->nreported,
			    th->hit[h]->nincluded,
			    (th->hit[h]->desc == NULL ? "-" : th->hit[h]->desc)) < 0)
                  ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-sequence hit list: write failed");
        }
    }
  return eslOK;
}


/* Function:  p7_tophits_TabularDomains()
 * Synopsis:  Output parseable table of per-domain hits
 *
 * Purpose:   Output a parseable table of reportable per-domain hits
 *            in sorted tophits list <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>, using final pipeline
 *            accounting stored in <pli>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEWRITE> if a write to <ofp> fails; for example, if
 *            the disk fills up.
 */
int
p7_etophits_TabularDomains(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, P7_PIPELINE *pli, int show_header)
{

  int qnamew = ESL_MAX(20, strlen(qname));
  int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int qaccw  = (qacc ? ESL_MAX(10, strlen(qacc)) : 10);
  int taccw  = ESL_MAX(10, p7_tophits_GetMaxAccessionLength(th));
  int tlen, qlen;
  int h,d,nd;

  if (show_header)
    {
      if (fprintf(ofp, "#%*s %30s %40s %11s %11s %11s\n", tnamew+qnamew-1+15+taccw+qaccw, "",                                   "-------- full sequence -------",        "-------------- this domain -------------",                "hmm coord",      "ali coord",     "env coord") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-domain hit list: write failed");
      if (fprintf(ofp, "#%-*s %-*s %5s %-*s %-*s %5s %9s %6s %5s %5s %5s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n",
		  tnamew-1, " target name",        taccw, "accession",  "tlen",  qnamew, "query name",           qaccw, "accession",  "qlen",  "E-value",   "score",  "bias",  "fwdsc", "spfwdsc", "time", "#",   "of",  "c-Evalue",  "i-Evalue",  "score",  "bias",  "from",  "to",    "from",  "to",   "from",   "to",    "acc",  "description of target") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-domain hit list: write failed");
      if (fprintf(ofp, "#%*s %*s %5s %*s %*s %5s %9s %6s %5s %5s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %4s %s\n", 
      tnamew-1, "-------------------", taccw, "----------", "-----", qnamew, "--------------------", qaccw, "----------", "-----", "---------", "------", "-----", "-----", "-----", "-----", "---", "---", "---------", "---------", "-----", "-----", "-----", "-----", "-----", "-----", "-----", "-----", "----", "---------------------") < 0)
        ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-domain hit list: write failed");
    }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
    {
        nd = 0;
        for (d = 0; d < th->hit[h]->ndom; d++)
          if (th->hit[h]->dcl[d].is_reported)
          {
              nd++;

              /* in hmmsearch, targets are seqs and queries are HMMs;
               * in hmmscan, the reverse.  but in the ALIDISPLAY
               * structure, lengths L and M are for seq and HMMs, not
               * for query and target, so sort it out.
               */
              if (pli->mode == p7_SEARCH_SEQS) { qlen = th->hit[h]->dcl[d].ad->M; tlen = th->hit[h]->dcl[d].ad->L;  }
              else                             { qlen = th->hit[h]->dcl[d].ad->L; tlen = th->hit[h]->dcl[d].ad->M;  }

              if (fprintf(ofp, "%-*s %-*s %5d %-*s %-*s %5d %9.2g %5.1f %6.1f %5.1f %5.1f %5.3f %3d %3d %9.2g %9.2g %6.1f %5.1f %5d %5d %5" PRId64 " %5" PRId64 " %5d %5d %4.2f %s\n",
			  tnamew, th->hit[h]->name,
			  taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
			  tlen,
			  qnamew, qname,
			  qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
			  qlen,
			  exp(th->hit[h]->lnP) * pli->Z,
			  th->hit[h]->score,
			  th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
			  th->hit[h]->fwd_score,
			  th->hit[h]->spfwd_score,
			  th->hit[h]->time,
			  nd,
			  th->hit[h]->nreported,
			  exp(th->hit[h]->dcl[d].lnP) * pli->domZ,
			  exp(th->hit[h]->dcl[d].lnP) * pli->Z,
			  th->hit[h]->dcl[d].bitscore,
			  th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* NATS to BITS at last moment */
			  th->hit[h]->dcl[d].ad->hmmfrom,
			  th->hit[h]->dcl[d].ad->hmmto,
			  th->hit[h]->dcl[d].ad->sqfrom,
			  th->hit[h]->dcl[d].ad->sqto,
			  th->hit[h]->dcl[d].iae,
			  th->hit[h]->dcl[d].ibe,
			  (th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].ibe - th->hit[h]->dcl[d].iae)))),
			  (th->hit[h]->desc ?  th->hit[h]->desc : "-")) < 0)
		ESL_EXCEPTION_SYS(eslEWRITE, "tabular per-domain hit list: write failed");
          }
    }
  return eslOK;
}


