/* ribosum_matrix.c
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

struct ribomatrix_s {
  ESL_DMATRIX *riboP;      // ribosum joint       16x16 matrix 
  ESL_DMATRIX *riboC;      // ribosum conditional 16x16 matrix
  ESL_DMATRIX *riboQ;      // ribosum rate        16x16 matrix
  ESL_DMATRIX *rnaP;       // unpaired joint        4x4 matrix
  ESL_DMATRIX *rnaC;       // unpaired conditional  4x4 matrix
  ESL_DMATRIX *rnaQ;       // unpaired rate         4x4 matrix
  
  double      *f;          /* background frequencies */
  
  char        *name;
}
  
extern int                  Ribosum_matrix_AddCounts(ESL_MSA *msa, struct ribomatrix_s *ribosum, float *thresh2);
extern int                  Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, char *errbuf);
extern struct ribomatrix_s *Ribosum_matrix_Create();
extern void                 Ribosum_matrix_Destroy(struct ribomatrix_s *ribosum);

#endif
