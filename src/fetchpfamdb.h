/* fetchpfamdb
 *
 */
#ifndef FETCHPFAMDB_INCLUDED
#define FETCHPFAMDB_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_tree.h"

#include "orthologs.h"

#define NTAXALEVELS 7	/* number of taxonomic leves used    */
enum taxotype_e {
  SUPERKINGDOM = 0,
  KINGDOM      = 1,
  PHYLUM       = 2,
  CLASS        = 3,
  ORDER        = 4,
  FAMILY       = 5,
  GENUS        = 6,
  UNKNOWN      = 7,
};

typedef struct {
  int    n;		/* number of sequences */

  char  **seqid;
  char  **spid;
  char  **spname;
  char ***parent;

  enum taxotype_e parentype;

  int   nalloc;	/* current allocated # of taxa */
  
} SPECIES;

enum myatype_e {
  DOMYA     = 0, /* calculate time in mya for all comparisons, store in distfile */
  FROMFILE  = 1, /* distfile exist */
  ONTHESPOT = 2,  /* calculate time in mya on the spot for each comparison */
};


extern SPECIES *species_Create(int n);
extern void     species_Destroy(SPECIES *SPE);
extern int      fetch_MSATreeAndSpecies(char **ret_msafile, char **ret_treefile, char **ret_spfile, char **ret_sptreefile, char **ret_distfile,  
					char *entry, int full, int group, enum taxotype_e taxotype, enum myatype_e myatype, char *errbuf);
extern int      fetch_ReadSpeciesFile(FILE *outfp, char *spfile, char *errbuf, SPECIES **ret_SPE, ESL_MSA *msa, int verbose);

#endif /* FETCHPFAMDB_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/
