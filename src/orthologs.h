/* orthologs
 *
 */
#ifndef ORTHOLGS_INCLUDED
#define ORTHOLGS_INCLUDED


#include <stdio.h>		/* FILE */

#include "easel.h"
#include "esl_msa.h"
#include "esl_random.h"

typedef struct {
  char **taxalist;         /* list of sequence names [0,..,ntaxa-1] */  
  int    ntaxa;            /* number of sequences */
  int    no;	 	   /* number of orthologs */

  double orthodist;	   /* orthology 'distance' (-log[rioP]) */

  int   *useme;            /* TRUE if ortholog [0,..,ns-1] */  
} ORTHO;

extern ORTHO *ortho_Create(int nc, char **taxalist, int ntaxa);
extern void   ortho_Destroy(int nc, ORTHO *ortho);
extern int    ortho_FindClusters(char **sptreefile, FILE *outfp, ESL_RANDOMNESS *r, ESL_MSA *msa, char *outheader, ORTHO **ret_ortho, int *ret_nc, int nboot, int reboot,
				 double minorthsc, char *errbuf, int verbose);
extern int    ortho_BoostrapTrees(ESL_RANDOMNESS *r, ESL_MSA *msa, char *outheader, int Nb, char *treebootfile, char *errbuf, int verbose);
extern char  *ortho_RunSEQBOOT(ESL_RANDOMNESS *r, ESL_MSA *msa, char *outheader, int N, char *errbuf, int verbose);
extern int    ortho_RunRIO(char **sptreefile, FILE *outfp, ESL_MSA *msa, char *treebootfile, ORTHO **ret_ortho, int *ret_nc, double minorthsc, char *errbuf, int verbose);
extern int    ortho_LinkageRIOTable(char *riotbl, FILE *outfp, ESL_MSA *msa, ORTHO **ret_ortho, int *ret_nc, double minorthsc, char *errbuf, int verbose);
extern int    ortho_CompleteLinkage(char *ssiriotbl, FILE *outfp, ESL_MSA *msa, char **riotaxalist, int riontaxa, ORTHO **ret_cortho, int *ret_nc, 
				    double minorthsc, char *errbuf, int verbose);
extern int    ortho_SingleLinkage(char *ssiriotbl, FILE *outfp, ESL_MSA *msa, char **riotaxalist, int riontaxa, ORTHO **ret_cortho, int *ret_nc, 
				  double minorthsc, char *errbuf, int verbose);
extern int    ortho_LinkageTree2Clusters(FILE *outfp, ESL_MSA *msa, ESL_TREE *T, ORTHO **ret_cortho, int *ret_nc, double minorthsc, char *errbuf,  int verbose);
#endif /* ORTHOLGS_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/
