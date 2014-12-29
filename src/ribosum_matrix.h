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

#define IDX(i,j,L)     ( (i) * (L) + (j) )

struct ribomatrix_s {
  ESL_DMATRIX *prnaP;      // paired joint         16x16 matrix: P(aa',bb') = P(bb',aa')
  ESL_DMATRIX *prnaC;      // paired conditional   16x16 matrix: P(aa'|bb') = P(aa',bb')/pm(bb')
  ESL_DMATRIX *prnaQ;      // paired rate          16x16 matrix: Q(aa',bb') = log P(aa'|bb')
  ESL_DMATRIX *urnaP;      // unpaired joint         4x4 matrix: P(a,b) = P(b,a)
  ESL_DMATRIX *urnaC;      // unpaired conditional   4x4 matrix: P(a|b) = P(a,b)/pm(b)
  ESL_DMATRIX *urnaQ;      // unpaired rate          4x4 matrix: Q(a,b) = log P(a|b)

  double      *prna;       // [16] marginals for paired   positions pm(aa') = \sum_{bb'}P(aa',bb')
  double      *urna;       // [4]  marginals for unpaired positions pm(a)   = \sum_{b}P(a,b)
  double      *bg;         // background frequencies
  
  ESL_ALPHABET *abc;
  char         *name;
};
  
extern int                  Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, FILE *fp,
						     double tol, int verbose, char *errbuf);
extern int                  Ribosum_matrix_ConditionalsFromJoint(struct ribomatrix_s *ribosum, double tol, int verbose, char *errbuf);
extern struct ribomatrix_s *Ribosum_matrix_Create(ESL_ALPHABET *abc, char *name);
extern void                 Ribosum_matrix_Destroy(struct ribomatrix_s *ribosum);
extern int                  Ribosum_matrix_RateFromConditionals(struct ribomatrix_s *ribosum, double tol, int verbose, char *errbuf);
extern int                  Ribosum_matrix_Write(FILE *fp, struct ribomatrix_s *ribosum);

#endif
