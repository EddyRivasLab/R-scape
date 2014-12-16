/*  fetchfamdb - wrapper around fetchPfamTreeAndSpecies.pl script
 * 
 * Contents:
 *   1. Miscellaneous functions for msatree
 *   2. Unit tests
 *   3. Test driver
 *   4. License and copyright 
 *
 * ER, Fri Oct 21 10:26:52 EDT 2011 [Janelia] 
 * SVN $Id:$
 */

#include "p7_config.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include <math.h>
#include <float.h>
	
#include "easel.h"
#include "esl_msa.h"
#include "esl_tree.h"

#include "hmmer.h"

#include "fetchpfamdb.h"
#include "ssifile.h"

/*****************************************************************
 * 1. functions for structure SPECIES
 *****************************************************************/
/* Function:  species_Create()
 *
 * Purpose:   Allocate an empty tree structure for <ntaxa> taxa
 *            and return a pointer to it. <ntaxa> must be $\geq 2$.
 *
 * Args:      <ntaxa>   - number of taxa
 *
 * Returns:   pointer to the new <ESL_TREE> object; caller frees 
 *            this with <esl_tree_Destroy()>.
 *
 * Throws:    <NULL> if allocation fails.
 */
SPECIES *
species_Create(int n)
{
  SPECIES *SPE = NULL;
  int       i;
  int       t;
  int       status;

  /* Contract verification  */
  ESL_DASSERT1((n >= 1));

  /* 1st allocation round  */
  ESL_ALLOC(SPE, sizeof(SPECIES));
  SPE->seqid  = NULL;
  SPE->spid   = NULL;
  SPE->spname = NULL;
  SPE->parent = NULL;
  
  /* 2nd allocation round */
  SPE->n = n;  
  SPE->parentype = UNKNOWN;
  ESL_ALLOC(SPE->seqid,       sizeof(char  *) * n);
  ESL_ALLOC(SPE->spid,        sizeof(char  *) * n);
  ESL_ALLOC(SPE->spname,      sizeof(char  *) * n);
  ESL_ALLOC(SPE->parent,      sizeof(char **) * n * NTAXALEVELS);
  for (i = 0; i < SPE->n; i++)
    ESL_ALLOC(SPE->parent[i], sizeof(char  *) * NTAXALEVELS);

  for (i = 0; i < SPE->n; i++)
    {
      SPE->seqid[i]  = NULL;
      SPE->spid[i]   = NULL;
      SPE->spname[i] = NULL;

      for (t = 0; t < NTAXALEVELS; t++)
	SPE->parent[i][t] = NULL;
   }

   return SPE;
  
 ERROR:
  species_Destroy(SPE);
  return NULL;
}

/* Function:  species_Destroy()
 *
 * Purpose:   Frees an <SPECIES> object.
 */
void
species_Destroy(SPECIES *SPE)
{
  int i;
  int t;

  if (SPE == NULL) return;

  if (SPE->seqid) {
    for (i = 0; i < SPE->n; i++) if (SPE->seqid[i]) free(SPE->seqid[i]);
    free(SPE->seqid);
  }
  if (SPE->spid) {
    for (i = 0; i < SPE->n; i++) if (SPE->spid[i]) free(SPE->spid[i]);
    free(SPE->spid);
  }
  if (SPE->spname) {
    for (i = 0; i < SPE->n; i++) if (SPE->spname[i]) free(SPE->spname[i]);
    free(SPE->spname);
  }
  if (SPE->parent) {
    for (i = 0; i < SPE->n; i++) {
      if (SPE->parent[i]) {
	for (t = 0; t < NTAXALEVELS; t ++)      
	  if (SPE->parent[i][t]) free(SPE->parent[i][t]);
	free(SPE->parent[i]); 
      }
    } 
    free(SPE->parent);
  }
  
  free(SPE);
  return;
}


/*****************************************************************
 * 2. Miscellaneous functions for fetchpfamdb
 *****************************************************************/
int
fetch_MSATreeAndSpecies(char **ret_msafile, char **ret_treefile, char **ret_spfile, char **ret_sptreefile, char **ret_distfile, char *entry, 
			int full, int group, enum taxotype_e taxotype, enum myatype_e myatype, char *errbuf)
{
  char *msafile    = NULL;
  char *treefile   = NULL;
  char *spfile     = NULL;
  char *sptreefile = NULL;
  char *distfile   = NULL;
  char *args       = NULL;
  char *optgroup   = NULL;
  char *s;
  char *st;
  char *acc = NULL;
  int   status;
  
  if (entry == NULL) { ESL_XFAIL(eslFAIL, errbuf, "cannot fin the acc for this alignment\n"); goto ERROR; }

  if ("EVOHMMDIR" == NULL)               { ESL_XFAIL(eslENOTFOUND, errbuf, "cannot fin env EVOHMDIR\n"); goto ERROR; }
  if ((s = getenv("EVOHMMDIR")) == NULL) { ESL_XFAIL(eslENOTFOUND, errbuf, "cannot fin env EVOHMDIR\n"); goto ERROR; }

  /* sometimes Pfam accesions have a version number after a dot .
   * For instance PF08702, has AC PF08702.5 */
  st = entry;
  esl_strtok(&st, ".", &acc);
   if (full) {
    esl_sprintf(&msafile,    "%s.full",   acc);
    esl_sprintf(&args,       "%s/scripts/fetchPfamTreeAndSpecies.pl  -entry %s -align full ", s, entry);
  }
  else {
    esl_sprintf(&msafile,    "%s.seed",   acc);
    esl_sprintf(&args,       "%s/scripts/fetchPfamTreeAndSpecies.pl -entry %s ", s, entry);
  }
  esl_sprintf(&treefile,   "%s.tree",     msafile);
  esl_sprintf(&spfile,     "%s.species",  msafile);
  esl_sprintf(&sptreefile, "%s.tree.xml", msafile);
  
  if (myatype == DOMYA) esl_strcat(&args, -1, "-time ", -1);
  switch(taxotype) {
  case UNKNOWN: 
    if (group > 0) {
      esl_sprintf(&optgroup, "-group %d ", group);
      esl_strcat(&args, -1, optgroup, -1);
      if (myatype != ONTHESPOT) esl_sprintf(&distfile,   "%s.group%d.distance", msafile, group);
    }
    free(optgroup);
    break;
  case SUPERKINGDOM: esl_strcat(&args, -1, "-level superkingdom ", -1); if (myatype != ONTHESPOT) esl_sprintf(&distfile,   "%s.superkingdom.distance", msafile);
    break;
  case KINGDOM:      esl_strcat(&args, -1, "-level kingdom ", -1);      if (myatype != ONTHESPOT) esl_sprintf(&distfile,   "%s.kingdom.distance", msafile);
    break;
  case PHYLUM:       esl_strcat(&args, -1, "-level phylum ", -1);       if (myatype != ONTHESPOT) esl_sprintf(&distfile,   "%s.phylum.distance", msafile);
    break;
  case CLASS:        esl_strcat(&args, -1, "-level class ", -1);        if (myatype != ONTHESPOT) esl_sprintf(&distfile,   "%s.class.distance", msafile);
    break;
  case ORDER:        esl_strcat(&args, -1, "-level order ", -1);        if (myatype != ONTHESPOT) esl_sprintf(&distfile,   "%s.order.distance", msafile);
    break;
  case FAMILY:       esl_strcat(&args, -1, "-level family ", -1);       if (myatype != ONTHESPOT) esl_sprintf(&distfile,   "%s.family.distance", msafile);
    break;
  case GENUS:        esl_strcat(&args, -1, "-level genus ", -1);        if (myatype != ONTHESPOT) esl_sprintf(&distfile,   "%s.genus.distance", msafile);
    break;
  default: ESL_XFAIL(eslFAIL, errbuf, "unknown taxonomy type \n"); goto ERROR;
  }

  /* if any of the files is missing, run script fetchPfamTreeAndSpecies.pl  */
  if (   !esl_FileExists(treefile)                          ||  
         !esl_FileExists(spfile)                            || 
         !esl_FileExists(sptreefile)                        || 
       ( !esl_FileExists(distfile) && myatype != ONTHESPOT)   ) system(args);
    
  if (ret_msafile)    *ret_msafile    = msafile;    else { remove(msafile);    free(msafile);    }
  if (ret_treefile)   *ret_treefile   = treefile;   else { remove(treefile);   free(treefile);   }
  if (ret_spfile)     *ret_spfile     = spfile;     else { remove(spfile);     free(spfile);     }
  if (ret_sptreefile) *ret_sptreefile = sptreefile; else { remove(sptreefile); free(sptreefile); }
  if (ret_distfile)   *ret_distfile   = distfile;   else {                     free(distfile);   }
  
  free(args);
  return eslOK;
  
 ERROR:
  if (args != NULL) free(args);
  remove(msafile);    if (msafile)    free(msafile); 
  remove(treefile);   if (treefile)   free(treefile);
  remove(spfile);     if (spfile)     free(spfile);
  remove(sptreefile); if (sptreefile) free(sptreefile);
                      if (distfile)   free(distfile);
  return status;  
}

int
fetch_ReadSpeciesFile(FILE *outf, char *spfile, char *errbuf, SPECIES **ret_SPE, ESL_MSA *msa, int verbose)
{
  char      *ssispfile = NULL;
  FILE      *spfp = NULL;
  SPECIES   *SPE = NULL;
  int        nseqt = 0; 
  int        nseq = 0;
  int        n;
  int        t;
  int        status;
 
  if (spfile == NULL) {
    SPE = species_Create(msa->nseq); 
    SPE->n = 0;
    *ret_SPE = SPE;
    return eslOK;
  }

  /* Create an ssi index of spfile.
   * spfile contains the info to associate species to Pfam sq names */
  status = ssi_spfile_Create(spfile, &ssispfile);
  if (status != eslOK) { printf("ssisp_Create failed for file %s\n", spfile); status = eslFAIL; goto ERROR; }
 
  if ((spfp = fopen(spfile, "r")) == NULL) { printf("failed to open %s\n", spfile); status = eslFAIL; goto ERROR;  }
  
  /* Create the tree, initially allocated for 2 taxa.
   * Allocate for taxon and node labels, too.
   */
  SPE = species_Create(msa->nseq); 
  if (SPE == NULL) { status = eslFAIL; goto ERROR; }
   
  for (n = 0; n < msa->nseq; n ++) {
    nseqt ++;
    esl_strdup(msa->sqname[n], -1, &(SPE->seqid[n]));
    ssi_spfile_GetSpeciesInfo(ssispfile, msa->sqname[n], n, SPE);  

    if (verbose) {
      printf(" taxon[%d] = %s | %s | %s ", n, SPE->seqid[n], SPE->spid[n], SPE->spname[n]);
      for (t = 0; t < NTAXALEVELS; t ++)
	printf("| %s ",  SPE->parent[n][t]);
      printf("\n");
    }
    if (SPE->spid[n] == NULL || SPE->spname[n]  == NULL || SPE->parent[n][0] == NULL)
      { fprintf(outf, "WARNING: couldn't find species match for sequence: %s\n", msa->sqname[n]); }
    else nseq ++;
  }
  SPE->n = msa->nseq;
  if (verbose) printf("PFAM/NCBI species match for %d / %d sequences\n\n", nseq, nseqt);

  fclose(spfp); 
  remove(ssispfile); free(ssispfile);
  *ret_SPE = SPE;
  return eslOK;

 ERROR:
  remove(ssispfile); if (ssispfile) free(ssispfile);
  if (spfp) fclose(spfp);
  if (SPE) species_Destroy(SPE);
  return status;

}
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
