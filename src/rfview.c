/* rfview.c */
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "rscape_config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "easel.h"

#include "correlators.h"
#include "rfview.h"

static int dir_exists(const char *path);

int 
rfview_Depict(char *rfviewfile, char *omsafile, char *covfile, int nagg, enum agg_e *agg_method, char *helixcovfile, int makepdf, int makesvg, char *errbuf, int verbose)
{
  char *rfviewsvg = NULL;
  char *rfviewpdf = NULL;
  char *package   = NULL;
  char *args      = NULL;
  int   status;

  if (!RSCAPE_BIN) ESL_XFAIL(status, errbuf, "Failed to find RFview executable\n");

  // Find the executable for the OS
  // MAC
#ifdef OS_MAC
  esl_sprintf(&package, "%s/../lib/RFview/RFview-%s/mac/RFview.app/Contents/MacOS/RFview", RSCAPE_BIN, RFVIEW_VERSION,
	      RSCAPE_BIN);
  
  // LINUX_X86_64
#elif defined(OS_LINUX_X86_64)
  esl_sprintf(&package, "%s/../src/squashfs-root/rfview ", RSCAPE_BIN);
  
  // LINUX_ARM64
#elif defined(OS_LINUX_ARM64)
  esl_sprintf(&package, "%s/../src/squashfs-root/rfview ", RSCAPE_BIN);

  // WIN - not supported
#elif defined(OS_WIN)
  printf("RFview for OS_WIN not implemented\n");
  
#else
  printf("OS not recognized\n");
#endif

  // so far only added MacOS executables
  if (!package) return eslOK;
  
  esl_sprintf(&rfviewpdf, "%s.pdf", rfviewfile);
  esl_sprintf(&rfviewsvg, "%s.svg", rfviewfile);

 
  if (nagg > 1) {
    printf("%s includes more than one aggregation method. RFview will not use it. Run each aggregation method separately.\n", helixcovfile);
    if (verbose) esl_sprintf(&args, "%s --incSsEnds --structureFile %s --basePairAnno %s --svg %s --pdf %s",
			     package, omsafile, covfile, rfviewsvg, rfviewpdf);
    else         esl_sprintf(&args, "%s --incSsEnds --structureFile %s --basePairAnno %s --svg %s --pdf %s >/dev/null",
			     package, omsafile, covfile, rfviewsvg, rfviewpdf);
  }
  else {
    if (agg_method[0] == AGG_NONE)
      if (verbose) esl_sprintf(&args, "%s --incSsEnds --structureFile %s --basePairAnno %s --svg %s --pdf %s",
			       package, omsafile, covfile, rfviewsvg, rfviewpdf);
      else         esl_sprintf(&args, "%s --incSsEnds --structureFile %s --basePairAnno %s --svg %s --pdf %s >/dev/null",
			       package, omsafile, covfile, rfviewsvg, rfviewpdf);
    else 
      if (verbose) esl_sprintf(&args, "%s --incSsEnds --structureFile %s --basePairAnno %s --helixCovAnno %s --svg %s --pdf %s",
			       package, omsafile, covfile, helixcovfile, rfviewsvg, rfviewpdf);
      else         esl_sprintf(&args, "%s --incSsEnds --structureFile %s --basePairAnno %s --helixCovAnno %s --svg %s --pdf %s >/dev/null",
			       package, omsafile, covfile, helixcovfile, rfviewsvg, rfviewpdf);
  }
  
 
  status = system(args);
  if (1||verbose) printf("%s\n", args);
  if (status == -1) ESL_XFAIL(status, errbuf, "Failed to run RFview\n");
  
  free(rfviewpdf);
  free(rfviewsvg);
  free(package);
  free(args);
  
  return eslOK;
  
 ERROR:
  if (rfviewpdf) free(rfviewpdf);
  if (rfviewsvg) free(rfviewsvg);
  if (package) free(package);
  if (args) free(args);
  
  return status;
}


static
int dir_exists(const char *path)
{
  struct stat stats;
  
  // stat() returns 0 on success
  if (stat(path, &stats) == 0) {
    // Check if the path is a directory
    return S_ISDIR(stats.st_mode);
  }
  
  return 0; // Does not exist or cannot be accessed
}
