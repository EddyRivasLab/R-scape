CFLAGS        = -O
CC            = cc

# ++++++++++++++++++++++++++++ Program RNADRAW
RNAVIEW = rnaview
OBJS_RNAVIEW = rnaview.o analyze.o fpair.o fpair_sub.o  pair_type.o  \
               multiple.o nrutil.o ps-xy.o  ps-xy-sub.o pattern.o  \
	       rnaxml-new.o statistics.o vrml.o xml2ps.o \

$(RNAVIEW):   $(OBJS_RNAVIEW)
	$(CC) -o $(RNAVIEW) $(OBJS_RNAVIEW) -lm

rnaview.o : rnaview.c
	cc -c  rnaview.c
fpair.o : fpair.c
	cc -c fpair.c
fpair_sub.o : fpair_sub.c
	cc -c fpair_sub.c
pair_type.o : pair_type.c
	cc -c pair_type.c
nrutil.o : nrutil.c 
	cc -c  nrutil.c
ps-xy.o  : ps-xy.c
	cc -c ps-xy.c
ps-xy-sub.o  : ps-xy-sub.c
	cc -c ps-xy-sub.c
vrml.o : vrml.c
	cc -c  vrml.c

rnaxml-new.o : rnaxml-new.c
	cc -c  rnaxml-new.c
analyze.o :  analyze.c
	cc -c  analyze.c
pattern.o :  pattern.c
	cc -c  pattern.c

xml2ps.o :  xml2ps.c
	cc -c  xml2ps.c

multiple.o :  multiple.c
	cc -c  multiple.c

statistics.o :  statistics.c
	cc -c  statistics.c


.PHONY : clean
clean:
	rm *.o
