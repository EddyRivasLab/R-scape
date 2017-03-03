/* ssifile
 * create a ssi file for a given species file
 * 
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_ssi.h"

#include "fetchpfamdb.h"
#include "ssifile.h"

int
ssi_spfile_Create(char *spfile, char **ret_ssispfile)
{
  char       *ssispfile = NULL;
  FILE       *spfp = NULL; 
  ESL_NEWSSI *ns = NULL;
  char       *buf = NULL;          /* growable buffer for esl_fgets()      */
  char       *seqid;
  char       *spid;
  char       *parent;
  char       *s;		
  uint16_t    fh;
  off_t       offset;
  int         n = 0;             /* length of buf                        */
  
  /* Create an ssi index of spfile */
  if (esl_strdup(spfile, -1, &ssispfile)     != eslOK) esl_fatal("ssi_Create() failed");
  if (esl_strcat(&ssispfile,  -1, ".ssi", 4) != eslOK) esl_fatal("ssi_Create() failed");
  if (esl_newssi_Open(ssispfile, TRUE, &ns)  != eslOK) esl_fatal("new SSI index open failed");
  
  if ((spfp = fopen(spfile, "r"))  == NULL) esl_fatal("failed to open %s\n", spfile); 
  if (esl_newssi_AddFile(ns, spfile, eslSQFILE_UNKNOWN, &fh) != eslOK) esl_fatal("esl_newssi_AddFile() failed");
  
  /* read through the file and add the keys */
  offset = ftello(spfp);
  while (esl_fgets(&buf, &n, spfp) == eslOK)
    {
      if (*buf == '#') { offset = ftello(spfp);  continue; } /* skip comments */

      s = buf;                      
      esl_strtok(&s, "\t", &spid);     /* species id    = 1st token online */
      esl_strtok(&s, "\t", &seqid);    /* sequence id = 2st token online */
      esl_strtok(&s, "\t", &parent);   /* parent        = 3st token online */
      if (esl_newssi_AddKey(ns, seqid, fh, offset, 0, 0) != eslOK)
	esl_fatal("failed to add key %s to index: %s", parent, ns->errbuf);
      offset = ftello(spfp);				 
    }
  free(buf);
  fclose(spfp);
  
  /* Save the SSI index to a file */
  esl_newssi_Write(ns);
  esl_newssi_Close(ns);
  
  *ret_ssispfile = ssispfile;
  return eslOK;
}

int 
ssi_spfile_GetSpeciesInfo(char *ssispfile, char *sqname, int n, SPECIES *SPE)
{
  ESL_SSI *ssi = NULL;
  uint16_t fh;                  /* file handle SSI associates w/ file   */
  char    *spfile = NULL;       /* name of spfile                       */
  int      fmt;                 /* format code (1, in this example)     */
  off_t    offset;              /* disk offset of seqname in spfile     */
  FILE    *fp = NULL;           /* opened spfile for reading            */
  char    *buf = NULL;          /* growable buffer for esl_fgets()      */
  char    *seqid = NULL;
  char    *spid = NULL;
  char    *spname = NULL;
  char    *parent = NULL;
  char    *sqnamedup = NULL;
  char    *seqidfetch = NULL;
  char    *s;
  int      t;
  int      size = 0;               /* size of buffer */
  int      status;

  if (ssispfile == NULL) return eslOK;
  if (sqname    == NULL) return eslOK;
  
  esl_strdup(sqname, -1, &sqnamedup);
  s = sqnamedup;
  esl_strtok(&s, "/", &seqid);  /* the seqid is what is before the / */
 
  if (esl_ssi_Open(ssispfile, &ssi)                          != eslOK) { printf("open failed\n"); status = eslFAIL; goto ERROR; }
  if (esl_ssi_FindName(ssi, seqid, &fh, &offset, NULL, NULL) != eslOK) {
    printf("\nWARNING: Could not find species for sqname %s\n", seqid); 
    status = eslOK;
    goto ERROR; /* exit with eslOK */
  }
  if (esl_ssi_FileInfo(ssi, fh, &spfile, &fmt)                  != eslOK) { printf("info failed\n"); status = eslFAIL; goto ERROR; }
  /* you can't close the ssi file yet - spfile is pointing into it! */

  if ((fp = fopen(spfile, "r"))     == NULL)  { printf("failed to open %s\n", spfile);     status = eslFAIL; goto ERROR; }
  if (fseeko(fp, offset, SEEK_SET)  != 0)     { printf("failed to position %s\n", spfile); status = eslFAIL; goto ERROR; }
  if (esl_fgets(&buf, &size, fp)       != eslOK) { printf("failed to get name/desc line\n");  status = eslFAIL; goto ERROR; }
  if (*buf != '#') { /* not a comment */
    s = buf;                      
    esl_strtok(&s, "\t", &spid);              /* spid       = 1st token online */
    esl_strtok(&s, "\t", &seqidfetch);        /* seqid      = 2st token online */
    for (t = 0; t < NTAXALEVELS; t++) {
      esl_strtok(&s, "!", &parent);          /* parent     = 3-9thth token online */
      esl_strdup(parent, -1, &SPE->parent[n][t]);
    }
    esl_strtok(&s, "\n", &spname);            /* spname     = remaining tokens online */  
    if (strcmp(seqidfetch, seqid) != 0) { printf("failed to find key %s\n", seqid); status = eslFAIL; goto ERROR; }
    
    esl_strdup(spid,   -1, &SPE->spid[n]);
    esl_strdup(spname, -1, &SPE->spname[n]);
  }
  fclose(fp);  
  esl_ssi_Close(ssi);  

  free(sqnamedup);
  free(buf);
  return eslOK;
  
 ERROR: 
  if (sqnamedup) free(sqnamedup);
  if (ssi) esl_ssi_Close(ssi);  
  if (buf) free(buf);
  return status;
}

int
ssi_distfile_Create(char *distfile, char **ret_ssidistfile, char ***ret_taxalist, int *ret_ntaxa)
{
  char       *ssidistfile = NULL;
  FILE       *distfp = NULL; 
  ESL_NEWSSI *ns = NULL;
  char       *buf = NULL;          /* growable buffer for esl_fgets()      */
  char       *species;
  char       *s;
  char      **taxalist = NULL;
  int         ntaxa = 0;
  uint16_t    fh;
  off_t       offset;
  int         n = 0;               /* length of buf                        */
  int         status;
  
  /* Create an ssi index of distfile */
  if (esl_strdup(distfile, -1, &ssidistfile)   != eslOK) esl_fatal("ssi_Create() failed");
  if (esl_strcat(&ssidistfile,  -1, ".ssi", 4) != eslOK) esl_fatal("ssi_Create() failed");
  if (esl_newssi_Open(ssidistfile, TRUE, &ns)  != eslOK) esl_fatal("new SSI index open failed");
  
  if ((distfp = fopen(distfile, "r"))  == NULL) esl_fatal("failed to open %s\n", distfile); 
  if (esl_newssi_AddFile(ns, distfile, eslSQFILE_UNKNOWN, &fh) != eslOK) esl_fatal("esl_newssi_AddFile() failed");
  
  /* allocate taxalist */
  ESL_ALLOC(taxalist, sizeof(char *) * (ntaxa+1)); 
  taxalist[0] = NULL;

  /* read through the file and add the keys */
  offset = ftello(distfp);
  while (esl_fgets(&buf, &n, distfp) == eslOK)
    {
      s = buf;                      
      esl_strtok(&s, "\t", &species);     /* species id    = 1st token online */

      if (esl_str_IsBlank(species)) { offset = ftello(distfp);	continue; } /* skip blanks */
    
      esl_strdealign(species, species, " ", NULL); /* remove possible blanks at the beginning */

      ESL_REALLOC(taxalist, sizeof(char *) * (ntaxa+1)); taxalist[ntaxa] = NULL;
      esl_strdup(species, -1, &(taxalist[ntaxa++]));
 
     if (esl_newssi_AddKey(ns, species, fh, offset, 0, 0) != eslOK)
	esl_fatal("failed to add key %s to index: %s", species, ns->errbuf);		
      offset = ftello(distfp);				 
    }
  free(buf);
  fclose(distfp);
	 
  /* Save the SSI index to a file */
  esl_newssi_Write(ns);
  esl_newssi_Close(ns);
  
  *ret_ssidistfile = ssidistfile;
  if (ret_taxalist) *ret_taxalist = taxalist; else free(taxalist);
  if (ret_ntaxa)    *ret_ntaxa    = ntaxa;
  return eslOK;

 ERROR:
  if (taxalist) free(taxalist);
  return status;
}

int 
ssi_distfile_GetDistance(char *ssidistfile, char **taxalist, int ntaxa, char *tx1, char *tx2, float *ret_mya)
{
  ESL_SSI  *ssi = NULL;
  uint16_t  fh;                  /* file handle SSI associates w/ file   */
  char     *distfile;            /* name of distfile                     */
  int       fmt;                 /* format code (1, in this example)     */
  off_t     offset;              /* disk offset of seqname in distfile   */
  FILE     *fp = NULL;           /* opened distfile for reading          */
  char     *buf = NULL;          /* growable buffer for esl_fgets()      */
  char     *s;
  char     *sp1;
  char     *sp2;
  char     *species = NULL;
  char     *distance = NULL;
  float     mya;
  int       x, x1, x2;
  int       X;
  int       n = 0;               /* size of buffer                       */
  int       status;

  if (ssidistfile == NULL) return eslOK;

  for (x = 0; x < ntaxa; x ++) {
     if (strcmp(tx1, taxalist[x]) == 0) { x1 = x; break; }
  }
  if (x == ntaxa) { printf("failed to find tx1 %s in taxalist\n", tx1); status = eslFAIL; goto ERROR; }
 
 for (x = 0; x < ntaxa; x ++) {
   if (strcmp(tx2, taxalist[x]) == 0) { x2 = x; break; }
  }
  if (x == ntaxa) { printf("failed to find tx2 %s in taxalist\n", tx2); status = eslFAIL; goto ERROR; }

  sp1 = (x1 < x2)? tx2 : tx1;
  sp2 = (x1 < x2)? tx1 : tx2;
  X   = (x1 < x2)? x1 : x2;
  
  if (esl_ssi_Open(ssidistfile, &ssi)                      != eslOK) { printf("open failed\n"); status = eslFAIL; goto ERROR; }
  if (esl_ssi_FindName(ssi, sp1, &fh, &offset, NULL, NULL) != eslOK) {
    printf("\nWARNING: Could not find taxon %s\n", sp1); return eslOK;
  }
  if (esl_ssi_FileInfo(ssi, fh, &distfile, &fmt)                  != eslOK) { printf("info failed\n"); status = eslFAIL; goto ERROR; }
  /* you can't close the ssi file yet - distfile is pointing into it! */

  if ((fp = fopen(distfile, "r"))   == NULL)  { printf("failed to open %s\n", distfile);     status = eslFAIL; goto ERROR; }
  if (fseeko(fp, offset, SEEK_SET)  != 0)     { printf("failed to position %s\n", distfile); status = eslFAIL; goto ERROR; }
  if (esl_fgets(&buf, &n, fp)       != eslOK) { printf("failed to get name/desc line\n");    status = eslFAIL; goto ERROR; }
  
  s = buf;                      
  esl_strtok(&s, "\t", &species);     /* species id    = 1st token online */

  /* make sure you got the right species */
  esl_strdealign(species, species, " ", NULL); /* remove possible blanks at the beginning */
  if (strcmp(species, sp1) != 0) {  printf("failed to find key %s but %s\n", sp1, species); status = eslFAIL; goto ERROR; }
  
  x = 0;
  while(x++ <= X) 
    esl_strtok(&s, "\t", &distance);  /* sp2 = X-th token online */
  
  esl_strdealign(distance, distance, " ", NULL); /* remove possible blanks at the beginning */
  if (strcmp(distance, "-") == 0) { mya = -1.0;           }
  else                            { mya = atof(distance); }   
  
  *ret_mya = mya;
  
  esl_ssi_Close(ssi);  
  fclose(fp);  
  free(buf);
  return eslOK;
  
 ERROR: 
  if (ssi) esl_ssi_Close(ssi);  
  if (fp) fclose(fp);
  if (buf) free(buf);
  return eslFAIL;
}

int
ssi_riotbl_Create(char *riotbl, char **ret_ssiriotbl, char ***ret_taxalist, int *ret_ntaxa)
{
  char       *ssiriotbl = NULL;
  FILE       *riofp = NULL; 
  ESL_NEWSSI *ns = NULL;
  char       *buf = NULL;          /* growable buffer for esl_fgets()      */
  char       *species;
  char       *boostrap;
  char       *s;
  char      **taxalist = NULL;
  int         ntaxa = 0;
  uint16_t    fh;
  off_t       offset;
  int         n = 0;               /* length of buf                        */
  int         status;
  
  /* Create an ssi index of distfile */
  if (esl_strdup(riotbl, -1, &ssiriotbl)     != eslOK) esl_fatal("ssi_Create() failed");
  if (esl_strcat(&ssiriotbl,  -1, ".ssi", 4) != eslOK) esl_fatal("ssi_Create() failed");
  if (esl_newssi_Open(ssiriotbl, TRUE, &ns)  != eslOK) esl_fatal("new SSI index open failed");
  
  if ((riofp = fopen(riotbl, "r"))  == NULL) esl_fatal("failed to open %s\n", riotbl); 
  if (esl_newssi_AddFile(ns, riotbl, eslSQFILE_UNKNOWN, &fh) != eslOK) esl_fatal("esl_newssi_AddFile() failed");
  
  /* allocate taxalist */
  ESL_ALLOC(taxalist, sizeof(char *) * (ntaxa+1)); 
  taxalist[0] = NULL;

  /* read through the file and add the keys */
  offset = ftello(riofp);
  while (esl_fgets(&buf, &n, riofp) == eslOK)
    {
      s = buf;                      
      esl_strtok(&s, "\t", &species);     /* species id    = 1st token online */
      esl_strtok(&s, "\t", &boostrap);    /* continues with all bootstrap values */

      if ( esl_str_IsBlank(species)) { offset = ftello(riofp);  continue; } /* skip blanks */

      if (!esl_str_IsReal(boostrap) && esl_strcmp(boostrap, "-") != 0) { offset = ftello(riofp);  continue; } /* ignore first line */

      esl_strdealign(species, species, " ", NULL); /* remove possible blanks at the beginning */
      ESL_REALLOC(taxalist, sizeof(char *) * (ntaxa+1)); taxalist[ntaxa] = NULL;
      esl_strdup(species, -1, &(taxalist[ntaxa++]));
  
      if (esl_newssi_AddKey(ns, species, fh, offset, 0, 0) != eslOK)
	esl_fatal("failed to add key %s to index: %s", species, ns->errbuf);		
      offset = ftello(riofp);				 
    }
  free(buf);
  fclose(riofp);
	 
#if 0
  for (n = 0; n < ntaxa; n ++) printf("taxa[%d]=%s\n", n, taxalist[n]);
#endif

  /* Save the SSI index to a file */
  esl_newssi_Write(ns);
  esl_newssi_Close(ns);
  
  *ret_ssiriotbl = ssiriotbl;
  if (ret_taxalist) *ret_taxalist = taxalist; else free(taxalist);
  if (ret_ntaxa)    *ret_ntaxa    = ntaxa;
  return eslOK;

 ERROR:
  if (taxalist) free(taxalist);
  return status;
}

int 
ssi_riotbl_GetOrthsc(char *ssiriotbl, char **taxalist, int ntaxa, char *tx1, char *tx2, float *ret_orthsc, int verbose)
{
  ESL_SSI  *ssi = NULL;
  uint16_t  fh;                  /* file handle SSI associates w/ file   */
  char     *riotbl;              /* name of riotbl                       */
  int       fmt;                 /* format code (1, in this example)     */
  off_t     offset;              /* disk offset of seqname in riotbl     */
  FILE     *fp = NULL;           /* opened riotbl for reading            */
  char     *buf = NULL;          /* growable buffer for esl_fgets()      */
  char     *s;
  char     *species = NULL;
  char     *boostrap = NULL;
  char     *sp1, *sp2;
  float     orthsc;
  int       x, x1, x2;
  int       X;
  int       n = 0;               /* size of buffer                       */
  int       status;

  if (ssiriotbl == NULL) return eslOK;

  /* it is possible those taxa were not analized by rio */
  for (x = 0; x < ntaxa; x ++) {
    if (strcmp(tx1, taxalist[x]) == 0) { x1 = x; break; }
  }
  if (x == ntaxa) { if (verbose) printf("RIO: does not analize %s\n", tx1); *ret_orthsc = 0.0; return eslOK; }
  
  for (x = 0; x < ntaxa; x ++) {
   if (strcmp(tx2, taxalist[x]) == 0) { x2 = x; break; }
  }
  if (x == ntaxa) { if (verbose) printf("RIO: does not analize %s\n", tx2); *ret_orthsc = 0.0; return eslOK; }

  sp1 = (x1 < x2)? tx2 : tx1;
  sp2 = (x1 < x2)? tx1 : tx2;
  X   = (x1 < x2)? x1  : x2;
  
  if (esl_ssi_Open(ssiriotbl, &ssi)                      != eslOK) { printf("open failed\n"); status = eslFAIL; goto ERROR; }
  if (esl_ssi_FindName(ssi, sp1, &fh, &offset, NULL, NULL) != eslOK) {
    printf("\nWARNING: Could not find taxon %s\n", sp1); return eslOK;
  }
  if (esl_ssi_FileInfo(ssi, fh, &riotbl, &fmt)                  != eslOK) { printf("info failed\n"); status = eslFAIL; goto ERROR; }
  /* you can't close the ssi file yet - riotbl is pointing into it! */

  if ((fp = fopen(riotbl, "r"))     == NULL)  { printf("failed to open %s\n", riotbl);     status = eslFAIL; goto ERROR; }
  if (fseeko(fp, offset, SEEK_SET)  != 0)     { printf("failed to position %s\n", riotbl); status = eslFAIL; goto ERROR; }
  if (esl_fgets(&buf, &n, fp)       != eslOK) { printf("failed to get name/desc line\n");  status = eslFAIL; goto ERROR; }
  
  s = buf;                      
  esl_strtok(&s, "\t", &species);     /* species id    = 1st token online */

  /* make sure you got the right species */
  esl_strdealign(species, species, " ", NULL); /* remove possible blanks at the beginning */
  if (strcmp(species, sp1) != 0) {  printf("failed to find key %s but %s\n", sp1, species); status = eslFAIL; goto ERROR; }
  
  x = 0;
  while(x++ <= X) 
    esl_strtok(&s, "\t", &boostrap);  /* sp2 = X-th token online */
  
  esl_strdealign(boostrap, boostrap, " ", NULL); /* remove possible blanks at the beginning */
  if (strcmp(boostrap, "-") == 0) { orthsc = 1.0;            }
  else                            { orthsc = atof(boostrap); }   

#if 0
  printf("%s / %s orthsc %f\n", tx1, tx2, orthsc);
#endif
  
  *ret_orthsc = orthsc;
  
  esl_ssi_Close(ssi);  
  fclose(fp);  
  free(buf);
  return eslOK;
  
 ERROR: 
  if (ssi) esl_ssi_Close(ssi);  
  if (fp) fclose(fp);
  if (buf) free(buf);
  return eslFAIL;
}

