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
#include "ratematrix.h"

typedef enum {
  STAR  = 0,         // independent sequences, star topology
  GIVEN = 1,         // use the tree determined by the input msa
  SIM   = 2          // simulate tree
} TREETYPE;


#endif /*COVSIMULATE_INCLUDED*/

/************************************************************
 * @LICENSE@
 ************************************************************/
