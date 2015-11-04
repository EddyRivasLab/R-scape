/* cov_simulate - funtions to generate synthetic msa with covariation
 *
 */
#ifndef COVSIMULATE_INCLUDED
#define COVSIMULATE_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_dmatrix.h"
#include "esl_histogram.h"
#include "esl_ratematrix.h"

#include "e2.h"
#include "e1_rate.h"
#include "ribosum_matrix.h"

#include "ratematrix.h"

typedef enum {
  STAR  = 0,         // independent sequences, star topology
  GIVEN = 1,         // use the tree determined by the input msa
  SIM   = 2          // simulate tree
} TREETYPE;

extern int cov_GenerateAlignment(ESL_RANDOMNESS *r, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
				 ESL_MSA **msafull, int noss, double tol, char *errbuf, int verbose);
extern int cov_GenerateAlignmentUngapped(ESL_RANDOMNESS *r, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
					 ESL_MSA **ret_msafull, int noss, double tol, char *errbuf, int verbose);
extern int cov_GenerateAlignmentAddGaps(ESL_RANDOMNESS *r, ESL_TREE *T, ESL_MSA *root, E1_RATE *e1rate, struct ribomatrix_s *ribosum, 
					ESL_MSA *msafull, int noss, double tol, char *errbuf, int verbose);

#endif /*COVSIMULATE_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
