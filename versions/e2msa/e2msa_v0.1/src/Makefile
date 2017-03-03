# VPATH and shell configuration
#
top_srcdir     = ..
srcdir         = .

SHELL          = /bin/sh


# Package information
#
PACKAGE         = e2msa
PACKAGE_VERSION = 0.1
PACKAGE_TARNAME = e2msa

E2MSA         = @E2MSA_NAME@
E2MSA_VERSION = 0.1.0
E2MSA_DATE    = MARCH 2017

# Installation targets
#
prefix      = /Users/rivase/src/src/mysource/versions/e2msa/e2msa_v0.1
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
bindir      = ${exec_prefix}/bin
libdir      = ${exec_prefix}/lib
includedir  = ${prefix}/include
mandir      = ${datarootdir}/man
docdir      = ${datarootdir}/doc/${PACKAGE_TARNAME}
pdfdir      = ${docdir}
mandir      = ${datarootdir}/man
man1dir     = ${mandir}/man1
man1ext     = .1

# Compiler configuration
#
CC        = gcc
CFLAGS    = -g -O2  
LDFLAGS   = -static 
SIMDFLAGS = @SIMD_CFLAGS@
CPPFLAGS  = 

# Other tools
#
AR        = /usr/bin/ar 
RANLIB    = ranlib
INSTALL   = /usr/bin/install -c


LIBS      = -leasel -lhmmer     -lm

IFLAGS   = -I../lib/hmmer/src -I../lib/hmmer/lib/easel 
LFLAGS   = -L../lib/hmmer/src -L../lib/hmmer/lib/easel 

# beautification magic stolen from git
#
QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 =
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_GEN     = @echo '    ' GEN $@;
	QUIET_AR      = @echo '    ' AR $@;
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR  $$subdir; \
		        ${MAKE} -C $$subdir
endif


# Code manifests
#
SUBDIRS = ../src

PROGS    = e1sim\
	   e2msa\
	   e2sim\
	   e2train\
           gappenalties\


OBJS     = e1_bg.o\
	   e1_model.o\
	   e1_rate.o\
	   e1_emit.o\
	   e1_goodness.o\
	   e1_simulate.o\
	   e2_evolve.o\
           e2_generic_decoding.o\
           e2_generic_fwdback.o\
           e2_generic_optacc.o\
           e2_generic_viterbi.o\
           e2_generic_vtrace.o\
           e2_gmx.o\
           e2_msa.o\
	   e2_pipeline.o\
           e2_profile.o\
           e2hmmer_profile.o\
           e2_profilesq.o\
           e2_trace.o\
           e2_train.o\
	   e2_tree.o\
           e2f_generic_decoding.o\
           e2f_generic_fwdback.o\
           e2f_generic_optacc.o\
           e2fhmmer_generic_decoding.o\
           e2fhmmer_generic_fwdback.o\
           e2fhmmer_generic_optacc.o\
           e2tracealign.o\
           evohmmer.o\
	   evotracealign.o\
           fastsp.o\
           fetchpfamdb.o\
           homology.o\
           logsum.o\
           minimize.o\
           miscellaneous.o\
           msatree.o\
           mya.o\
           msaprobs.o\
           msamanip.o\
           muscle.o\
	   modelconfig.o\
           ncbiblast.o\
           orthologs.o\
           phmmer.o\
           plot.o\
           p7_evopipeline.o\
           etophits_output_tabular.o\
           ratematrix.o\
           ratebuilder.o\
           fasta.o\
           ssifile.o\
	   tkf_model.o\
	   tkf_rate.o\
	   tkf_simulate.o\
	   ribosum_matrix.o\

all: ${PROGS} 


${PROGS}: %: %.o ${OBJS}
	${QUIET_CC}${CC} ${CFLAGS} ${IFLAGS} ${LFLAGS} -o $@ $@.o ${OBJS} ${LIBS}


.c.o:   
	${QUIET_CC}${CC} ${IFLAGS} ${CFLAGS} -o $@ -c $<

libsrc.a:  ${OBJS}
	${QUIET_AR}${AR} -r libsrc.a $(OBJS) > /dev/null 2>&1
	@${RANLIB} libsrc.a

install:
	${QUIET}if [ ! -z "${PROGS}" ]; then \
	   for file in ${PROGS}; do \
	      echo '    ' INSTALL $$file ;\
	      ${INSTALL} -m 0755 $$file ${DESTDIR}${bindir}/ ;\
	   done ;\
	fi


clean:
	-rm -f ${PROGS} ${OBJS} libsrc.a *.tab.c *.output *.o programs/*.o
	-rm -f *~ TAGS
	-rm -f *.gcno *.gcda *.gcov
	-rm -f core.[0-9]*
	-rm -f config.log config.status

distclean:
	-rm -f ${PROGS} ${OBJS} libsrc.a *.tab.c *.output *.o programs/*.o
	-rm -f *~ TAGS
	-rm -f *.gcno *.gcda *.gcov
	-rm -f core.[0-9]*
	-rm -f config.log config.status
	-rm -rf autom4te.cache
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out
	-rm -f cscope.po.out cscope.out cscope.in.out cscope.files
	-rm Makefile



################################################################
# @LICENSE@
################################################################
