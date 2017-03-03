/* Routines for the E2_PROFILE structure - e2 profile with two evolutionary models
 *                                         
 *    1. The E2_PROFILE object: allocation, initialization, destruction.
 *    2. Access methods.
 *    3. Debugging and development code.
 *    4. Unit tests.
 *    5. Test driver.
 *
 * See also: 
 *   modelconfig.c : routines that configure a profile given two evolutionary models
 */
#include "p7_config.h"
#include "hmmer.h"

#include "esl_vectorops.h"

#include "e2.h"
#include "e2hmmer_profile.h"


/* Function:  e2hmmer_profile_Create()
 * Synopsis:  Allocates a profile.
 *
 * Purpose:   Allocates for a profile of up to <M> nodes, for digital
 *            alphabet <abc>.
 *            
 *            Because this function might be in the critical path (in
 *            hmmscan, for example), we leave much of the model
 *            unintialized, including scores and length model
 *            probabilities. The <e2hmmer_ProfileConfig()> call is what
 *            sets these. 
 *            
 *            The alignment mode is set to <e2hmmer_NOMODE>.  The
 *            reference pointer <gm->abc> is set to <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation error.
 *
 * Xref:      STL11/125.
 */
E2HMMER_PROFILE *
e2hmmer_profile_Create(int M, const ESL_ALPHABET *abc)
{
  E2HMMER_PROFILE *gm = NULL;
  int              k;
  int              x;
  int              status;

  /* level 0 */
  ESL_ALLOC(gm, sizeof(E2HMMER_PROFILE));

  /* level 1 */
  ESL_ALLOC(gm->tsc,  sizeof(float  *) * (M+1) * e2HP_NTRANS);
  ESL_ALLOC(gm->sssc, sizeof(float **) * (M+1) * e2_MAXABET * e2_MAXABET);
  ESL_ALLOC(gm->ssc,  sizeof(float **) * (M+1) * e2_MAXABET * e2HP_NS);
  ESL_ALLOC(gm->isc,  sizeof(float **) * (M+1) * e2_MAXABET * e2HP_NS);
  
  for (k = 0; k <= M; k ++) {
    ESL_ALLOC(gm->tsc[k],  sizeof(float  ) * e2HP_NTRANS);
    ESL_ALLOC(gm->sssc[k], sizeof(float *) * e2_MAXABET * e2_MAXABET);
    ESL_ALLOC(gm->ssc[k],  sizeof(float *) * e2_MAXABET * e2HP_NS);
    ESL_ALLOC(gm->isc[k],  sizeof(float *) * e2_MAXABET * e2HP_NS);
    
    for (x = 0; x < abc->K;   x++) {
      ESL_ALLOC(gm->sssc[k][x], sizeof(float) * e2_MAXABET);
      ESL_ALLOC(gm->ssc[k][x],  sizeof(float) * e2HP_NS);
      ESL_ALLOC(gm->isc[k][x],  sizeof(float) * e2HP_NS);
    }
  }

  /* Set remaining info  */
  gm->mode = e2_NOMODE;
  gm->name = NULL;
  gm->acc  = NULL;
  gm->desc = NULL;  
  gm->abc  = abc;
  gm->M    = M;
  return gm;

 ERROR:
  e2hmmer_profile_Destroy(gm);
  return NULL;
}


/* Function:  e2hmmer_profile_Copy()
 * Synopsis:  Copy a profile.
 *
 * Purpose:   Copies profile <src> to profile <dst>, where <dst>
 *            has already been allocated to be of sufficient size.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <dst> is too small 
 *            to fit <src>.
 */
int
e2hmmer_profile_Copy(const E2HMMER_PROFILE *src, E2HMMER_PROFILE *dst)
{
  int k;
  int x;
  int status;

  dst->M = src->M;

  for (k = 0; k <= src->M; k ++) {
    esl_vec_FCopy(src->tsc[k], e2HP_NTRANS, dst->tsc[k]);
    for (x = 0; x < src->abc->K;   x++) esl_vec_FCopy(src->sssc[k][x], e2_MAXABET, dst->sssc[k][x]);
    for (x = 0; x < src->abc->K;   x++) esl_vec_FCopy(src->ssc[k][x],  e2HP_NS,     dst->ssc[k][x]);
    for (x = 0; x < src->abc->K;   x++) esl_vec_FCopy(src->isc[k][x],  e2HP_NS,     dst->isc[k][x]);
  }
  
  dst->mode = src->mode;

  if (dst->name != NULL) free(dst->name);
  if (dst->acc  != NULL) free(dst->acc);
  if (dst->desc != NULL) free(dst->desc);

  if ((status = esl_strdup(src->name,      -1, &(dst->name)))      != eslOK) return status;
  if ((status = esl_strdup(src->acc,       -1, &(dst->acc)))       != eslOK) return status;
  if ((status = esl_strdup(src->desc,      -1, &(dst->desc)))      != eslOK) return status;

  return eslOK;
}

/* Function:  e2hmmer_profile_Clone()
 * Synopsis:  Duplicates a profile.
 *
 * Purpose:   Duplicate profile <gm>; return a pointer
 *            to the newly allocated copy.
 */
E2HMMER_PROFILE *
e2hmmer_profile_Clone(const E2HMMER_PROFILE *gm)
{
  E2HMMER_PROFILE *g2 = NULL;
  int         status;

  if ((g2 = e2hmmer_profile_Create(gm->M, gm->abc)) == NULL) return NULL;
  if ((status = e2hmmer_profile_Copy(gm, g2)) != eslOK) goto ERROR;
  return g2;
  
 ERROR:
  e2hmmer_profile_Destroy(g2);
  return NULL;
}


/* Function:  e2hmmer_profile_Reuse()
 * Synopsis:  Prepare profile to be re-used for a new HMM.
 *
 * Purpose:   Prepare profile <gm>'s memory to be re-used
 *            for a new HMM.
 */
int
e2hmmer_profile_Reuse(E2HMMER_PROFILE *gm)
{
  /* name, acc, desc annotation is dynamically allocated for each HMM */
  if (gm->name != NULL) { free(gm->name); gm->name = NULL; }
  if (gm->acc  != NULL) { free(gm->acc);  gm->acc  = NULL; }
  if (gm->desc != NULL) { free(gm->desc); gm->desc = NULL; }

  /* reset some other things, but leave the rest alone. */
  gm->mode = e2_NOMODE;

  return eslOK;
}


/* Function:  e2hmmer_profile_Destroy()
 * Synopsis:  Frees a profile.
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 *
 * Xref:      STL11/125.
 */
void
e2hmmer_profile_Destroy(E2HMMER_PROFILE *gm)
{
  int k;
  int x;

  if (gm) {
    if (gm->tsc)  {
      for (k = 0; k <= gm->M; k ++) free(gm->tsc[k]);
      free(gm->tsc);
    }
    if (gm->sssc) {
      for (k = 0; k <= gm->M; k ++) {
	for (x = 0; x < gm->abc->K; x ++) free(gm->sssc[k][x]);
	free(gm->sssc[k]);
      }
      free(gm->sssc);
    }
    if (gm->ssc)  {
      for (k = 0; k <= gm->M; k ++) {
	for (x = 0; x < gm->abc->K; x ++) free(gm->ssc[k][x]);
	free(gm->ssc[k]);
      }
      free(gm->ssc);
    }
    if (gm->isc)  {
      for (k = 0; k <= gm->M; k ++) {
	for (x = 0; x < gm->abc->K; x ++) free(gm->isc[k][x]);
	free(gm->isc[k]);
      }
      free(gm->isc);
    }
    
    if (gm->name) free(gm->name);
    if (gm->acc)  free(gm->acc);
    if (gm->desc) free(gm->desc);
    free(gm);
  }
  return;
}

