/* ribosum_matrix.h
 *
 *   
*/
#ifndef RIBOSUM_MATRIX_INCLUDED
#define RIBOSUM_MATRIX_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"
#include "esl_alphabet.h"	/* ESL_DSQ, ESL_ALPHABET */
#include "esl_dmatrix.h"	/* ESL_DMATRIX           */

#define IDX(i,j,L)  ( (i) * (L) + (j) )

struct ribomatrix_s {
  char         *name;
  ESL_ALPHABET *abc;

  ESL_DMATRIX *bprsJ;      // basepairs joint                        16x16 matrix: bprsJ(aa',bb') = bprsJ(bb',aa')
  ESL_DMATRIX *bprsC;      // basepairs conditional                  16x16 matrix: bprsC(aa'|bb') = bprsJ(aa',bb')/bprsM(bb')
  ESL_DMATRIX *bprsQ;      // basepairs rate                         16x16 matrix: bprsQ(aa',bb') = log bprsC(aa'|bb')
  double      *bprsM;      // [16] marginals for basepairs                         bprsM(aa')     = \sum_{bb'} bprsJ(aa',bb')

  ESL_DMATRIX *prnaJ;      // paired positions joint                   4x4 matrix:     prnaJ(a,b) = prnaJ(b,a)
  ESL_DMATRIX *prnaC;      // paired positions conditional             4x4 matrix:     prnaC(a|b) = prnaJ(a,b)/prnaM(b)
  ESL_DMATRIX *prnaQ;      // paired positions rate                    4x4 matrix:     prnaQ(a,b) = log prnaC(a|b)
  double      *prnaM;      // [4]  marginals for paired positions                      prnaM(a)   = \sum_{b} prnaJ(a,b)

  ESL_DMATRIX *urnaJ;      // unpaired positions joint                 4x4 matrix:     urnaJ(a,b) = urnaJ(b,a)
  ESL_DMATRIX *urnaC;      // unpaired positions conditional           4x4 matrix:     urnaC(a|b) = urnaJ(a,b)/urnaM(b)
  ESL_DMATRIX *urnaQ;      // unpaired positions rate                  4x4 matrix:     urnaQ(a,b) = log urnaC(a|b)
  double      *urnaM;      // [4]  marginals for unpaired positions                    urnaM(a)   = \sum_{b} urnaJ(a,b)

  ESL_DMATRIX *xrnaJ;      // allpositions joint                       4x4 matrix:     xrnaJ(a,b) = xrnaJ(b,a)
  ESL_DMATRIX *xrnaC;      // allpositions conditional                 4x4 matrix:     xrnaC(a|b) = xrnaJ(a,b)/xrnaM(b)
  ESL_DMATRIX *xrnaQ;      // allpositions rate                        4x4 matrix:     xrnaQ(a,b) = log xrnaC(a|b)
  double      *xrnaM;      // [4]  marginals for allpositions                          xrnaM(a)   = \sum_{b} xrnaJ(a,b)

  double      *bg;         // background frequencies
};
  
extern int                  Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, FILE *fp,
						     double tol, int verbose, char *errbuf);;
extern int                  Ribosum_matrix_JointsFromMSA(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, 
							 double tol, int verbose, char *errbuf)
extern int                  Ribosum_matrix_ConditionalsFromJoint(struct ribomatrix_s *ribosum, double tol, int verbose, char *errbuf);
extern struct ribomatrix_s *Ribosum_matrix_Create(ESL_ALPHABET *abc, char *name);
extern void                 Ribosum_matrix_Destroy(struct ribomatrix_s *ribosum);
extern int                  Ribosum_matrix_RateFromConditionals(struct ribomatrix_s *ribosum, double tol, int verbose, char *errbuf);
extern int                  Ribosum_matrix_Write(FILE *fp, struct ribomatrix_s *ribosum);

#endif
