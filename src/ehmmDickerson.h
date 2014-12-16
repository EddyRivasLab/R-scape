#ifndef ehmmDICKERSON_INCLUDED
#define ehmmDICKERSON_INCLUDED

#include <stdio.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_fileparser.h>
#include <esl_msafile.h>
#include <esl_random.h>

#include "e2.h"
#include "p7_evopipeline.h"
#include "fetchpfamdb.h"
#include "plot.h"


/* struct cfg_s : "Global" application configuration shared by all threads/processes.
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  /* Shared configuration in masters & workers */
  int              argc;
  char           **argv;
  
  char             errbuf[eslERRBUFSIZE];
  ESL_RANDOMNESS  *r;	               /* random numbers for stochastic sampling (grm-emit) */
  ESL_ALPHABET    *abc;                /* the alphabet */
  double           fragfrac;	       /* seqs less than x*avg length are removed from alignment  */
  double           idthresh;	       /* fractional identity threshold for selecting subset of sequences */

  char            *msafile;
  char            *pfamname;
  char            *gnuplot;
  FILE            *outfp; 
  char            *outheader;          /* header for all output files */

  char            *hmmfile;            /* for e2hmmer case: open input HMM file */

  P7_PIPELINE     *pli;
  P7_BG           *bg;	   	       /* null model (copies made of this into threads) */
 
  P7_RATE         *R;

  int              nr;
  float            nr_scale;
  char            *subsmx;             /* BLUSUM62, ... */
  float            scaleR;             /* scale substitutions rate matrix. if 0, do not rescale */
 
  EVOM             evomodel;
  float            betainf;
  float            fixtime;

  char            *xlabel;
  char            *ylabel;

  int              nboot;               /* number of boostrapped genes trees */
  int              reboot;              /* refresh the treeboot file even if it exists */

  int              clsize;              /* minimum orthologs cluster size to analyze */
  int              clmax;               /* if true, analyze the maximal cluster of orthologs */
  float            orthcutoff;
  float            evalcutoff;

  int              submsa;              /* set to the number of random seqs taken from original msa.
					 * Set to 0 if we are taking all */
  int              infmt;
  int              noopt;               /* TRUE for no optimization of parameters */
  int              makeT;
  int              full;                /* TRUE if alignment is full (default seed) */
  int              group;               /* the approximate number of groups the seqs should be grouped into */
  enum taxotype_e  taxotype;            /* supercedes group, allows to specify a taxonomy level to group by,
			                 * options are (superkingdom, kingdom, phylum, class, order, family, genus) */

  int              do_viterbi;

  float            bintime;             /* time interval for binning data in binned plots */
  float            tlinear;             /* max time for a linear approximation. if 0, do all points */
  enum myatype_e   myatype;             /* how are we obtaining the time in mya for each species comparison */

  enum plottype_e  plottype;
  int              view;

  float            tol;
  int              verbose;
};


#endif
