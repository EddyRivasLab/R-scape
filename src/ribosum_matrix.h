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
  ESL_DMATRIX *prnaP;      // paired joint       16x16 matrix 
  ESL_DMATRIX *prnaC;      // paired conditional 16x16 matrix
  ESL_DMATRIX *prnaQ;      // paired rate        16x16 matrix
  ESL_DMATRIX *urnaP;      // unpaired joint        4x4 matrix
  ESL_DMATRIX *urnaC;      // unpaired conditional  4x4 matrix
  ESL_DMATRIX *urnaQ;      // unpaired rate         4x4 matrix
  double      *bg;         // background frequencies
  
  char        *name;
}
  
extern int                  Ribosum_matrix_Calculate(ESL_MSA *msa, struct ribomatrix_s *ribosum, float thresh1, float thresh2, char *errbuf);
extern struct ribomatrix_s *Ribosum_matrix_Create();
extern void                 Ribosum_matrix_Destroy(struct ribomatrix_s *ribosum);

#endif
