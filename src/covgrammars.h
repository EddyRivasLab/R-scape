/* covgrammarsr.h
 *
 *   
*/
#ifndef COVGRAMMARS_INCLUDED
#define COVGRAMMARS_INCLUDED

#include "p7_config.h"

#include <stdio.h>		/* FILE */


#include "easel.h"
#include "esl_random.h"


#define MAXLOOP_H 30   // maximum loop length
#define MAXLOOP_B 30   // maximum loop length
#define MAXLOOP_I 30   // maximum loop length
#define MINHAIRPIN 5   // minumum length of a hairping (including the 2 closing pairs)

enum grammar_e {
  G6,
  G6S,
  BGR
};

typedef float SCVAL;

typedef struct {
  SCVAL  t1[2];  // S -> LS     | L
  SCVAL  t2[2];  // L -> a F a' | a
  SCVAL  t3[2];  // F -> a F a' | LS

  SCVAL e_sing[4];
  SCVAL e_pair[16];
} G6param;

typedef struct {
  SCVAL  t1[2];  // S -> LS     | L
  SCVAL  t2[2];  // L -> a F a' | a
  SCVAL  t3[2];  // F -> a F a' | Ls

  SCVAL e_sing[4];
  SCVAL e_pair[16];
  SCVAL e_stck[16][16];
} G6Sparam;

typedef struct {
  SCVAL tS[3];     // S  -> a s | F0 S | epsilon
  SCVAL tF0[2];    // F0 -> a F5 a' | a P a'
  SCVAL tF5[2];    // F5 -> a F5 a' | a P a'
  SCVAL tP[5];     // P  -> m..m | m..m F0 | F0 m..m | d..d F0 d..d | M1 M
  SCVAL tM[2];     // M  -> M M1 | R
  SCVAL tR[2];     // R  ->  R a | M1
  SCVAL tM1[2];    // M1 -> a M1 | F0

  SCVAL e_sing[4];
  SCVAL e_pair1[16];
  SCVAL e_pair2[16];
  SCVAL e_stck1[16][16];
  SCVAL e_stck2[16][16];

  SCVAL e_single_l1[4];
  SCVAL l1[MAXLOOP_H]; // hairpin  loops
  
  SCVAL e_single_l2[4];
  SCVAL l2[MAXLOOP_B]; // bulge    loops
  
  SCVAL e_single_l3[4];
  SCVAL l3[MAXLOOP_I]; // internal loops
} BGRparam;



typedef struct {
  int       L;    /* msa length */
  SCVAL  **dp;   /* L * L triangular DP matrix */
} GMX;

typedef struct {
  GMX *S;
  GMX *L;
  GMX *F;
} G6_MX;

typedef struct {
  GMX *S;
  GMX *F0;
  GMX *F5;
  GMX *P;
  GMX *M;
  GMX *R;
  GMX *M1;

} BGR_MX;

extern GMX    *GMX_Create   (int L);
extern G6_MX  *G6MX_Create  (int L);
extern BGR_MX *BGRMX_Create (int L);
extern void    GMX_Destroy  (GMX *gmx);
extern void    G6MX_Destroy (G6_MX *g6mx);
extern void    BGRMX_Destroy(BGR_MX *bgrmx);
extern void    GMX_Dump     (FILE *fp, GMX *gmx);


#endif
