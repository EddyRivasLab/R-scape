/* src/rscape_config.h.  Generated from rscape_config.h.in by configure.  */
/* @configure_input@
 * rscapeconfig.h.in -> rscapeconfig.h
 * 
 * rscapeconfig.h is generated from rscapeconfig.h.in by the ./configure script.
 * DO NOT EDIT rscapeconfig.h; only edit rscapeconfig.h.in.
 *
 * ER, Thu Dec  8 11:53:34 EST 2011 [janelia] 
 * SVN $Id: rscape_config.h.in  $
 */
#ifndef RSCAPE_CONFIGH_INCLUDED
#define RSCAPE_CONFIGH_INCLUDED

/* The symbol alphabet is handled by ESL_ALPHABET objects, which
 * dynamically allocate; but sometimes I use statically-allocated
 * space, and it's useful to know a reasonable maximum for
 * symbol alphabet size.
 */
#define rscape_MAXABET    20      /* maximum size of alphabet (4 or 20)              */
#define rscape_MAXCODE    29      /* maximum degenerate alphabet size (18 or 29)     */

/* Version info - set once for whole package in configure.ac
 */
#define RSCAPE_VERSION "0.01"
#define RSCAPE_DATE "NOV 2015"
#define RSCAPE_COPYRIGHT "Copyright (C) 2015 Howard Hughes Medical Institute."
#define RSCAPE_LICENSE "Freely distributed under the GNU General Public License (GPLv3)."
#define RSCAPE_URL "http://hmmer.org/"

/* Large file support (must precede any header file inclusion.)
 */
/* #undef _FILE_OFFSET_BITS */
/* #undef _LARGE_FILES */
/* #undef _LARGEFILE_SOURCE */

/* Choice of optimized implementation (one and only one must be set)
 */
/* #undef rscape_IMPL_SSE */
/* #undef rscape_IMPL_VMX */
/* #undef rscape_IMPL_DUMMY */

/* Optional parallel implementations
 */
/* #undef HAVE_SSRSCAPE */
/* #undef HAVE_MPI */
/* #undef RSCAPE_PVM */
/* #undef EHMM_THREADS */
#define RSCAPE_THREADS 1
/* #undef HAVE_PTHREAD_ATTR_SETSCOPE */
/* #undef HAVE_PTHREAD_SETCONCURRENCY */

/* Optional processor specific support
 */
/* #undef HAVE_FLUSH_ZERO_MODE */

/* Debugging hooks
 */
/* #undef rscape_DEBUGGING */

#endif /*RSCAPE_CONFIGH_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
