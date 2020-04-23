CC        = cc



#----------------------------------------------------------------------------
# Project specific path defintions.
#----------------------------------------------------------------------------
PROJDIR    = .

SRC        = $(PROJDIR)/src
INCL       = $(PROJDIR)/include
OBJ        = $(PROJDIR)/obj
BIN        = $(PROJDIR)/bin
TEST       = $(PROJDIR)/test

LDEFINES    = -DNO_RANGE_CHECK
LINCLUDES      =  -I$(INCL)

#----------------------------------------------------------------------------
# Target files
#----------------------------------------------------------------------------
RNAVIEW    = $(BIN)/rnaview

SRCFILES = $(SRC)/rnaview.c \
           $(SRC)/fpair.c  \
           $(SRC)/fpair_sub.c  \
           $(SRC)/pair_type.c  \
           $(SRC)/nrutil.c  \
           $(SRC)/ps-xy.c  \
           $(SRC)/ps-xy-sub.c  \
           $(SRC)/vrml.c  \
           $(SRC)/rnaxml-new.c  \
           $(SRC)/analyze.c   \
           $(SRC)/pattern.c  \
           $(SRC)/xml2ps.c  \
           $(SRC)/multiple.c \
           $(SRC)/statistics.c


HFILES = $(INCL)/rna.h $(INCL)/nrutil.h $(INCL)/rna_header.h \
	$(INCL)/vrml.h $(INCL)/xml2ps.h


OBJ_FILE = $(OBJ)/rnaview.o \
           $(OBJ)/fpair.o  \
           $(OBJ)/fpair_sub.o  \
           $(OBJ)/pair_type.o  \
           $(OBJ)/nrutil.o  \
           $(OBJ)/ps-xy.o  \
           $(OBJ)/ps-xy-sub.o  \
           $(OBJ)/vrml.o  \
           $(OBJ)/rnaxml-new.o  \
           $(OBJ)/analyze.o   \
           $(OBJ)/pattern.o  \
           $(OBJ)/xml2ps.o  \
           $(OBJ)/multiple.o \
           $(OBJ)/statistics.o

all: $(RNAVIEW)

CFLAGS  =  $(LINCLUDES) 

$(RNAVIEW) : $(HFILES) $(OBJ_FILE) 
	$(CC) $(CFLAGS) -o $@ $(OBJ_FILE) $(LDFLAGS) -lm $(MALLOCLIB)



$(OBJ)/rnaview.o : $(SRC)/rnaview.c 
	$(CC) $(CFLAGS) -c $(SRC)/rnaview.c -o $@

$(OBJ)/fpair.o : src/fpair.c 
	$(CC) $(CFLAGS) -c $(SRC)/fpair.c -o $@

$(OBJ)/fpair_sub.o : $(SRC)/fpair_sub.c
	$(CC) $(CFLAGS) -c $(SRC)/fpair_sub.c -o $@

$(OBJ)/pair_type.o : $(SRC)/pair_type.c 
	$(CC) $(CFLAGS) -c $(SRC)/pair_type.c -o $@

$(OBJ)/nrutil.o : $(SRC)/nrutil.c 
	$(CC) $(CFLAGS) -c  $(SRC)/nrutil.c -o $@

$(OBJ)/ps-xy.o  : $(SRC)/ps-xy.c 
	$(CC) $(CFLAGS) -c $(SRC)/ps-xy.c -o $@

$(OBJ)/ps-xy-sub.o  : $(SRC)/ps-xy-sub.c 
	$(CC) $(CFLAGS) -c $(SRC)/ps-xy-sub.c -o $@

$(OBJ)/vrml.o : $(SRC)/vrml.c 
	$(CC) $(CFLAGS) -c  $(SRC)/vrml.c -o $@

$(OBJ)/rnaxml-new.o : $(SRC)/rnaxml-new.c 
	$(CC) $(CFLAGS) -c  $(SRC)/rnaxml-new.c -o $@

$(OBJ)/analyze.o :  $(SRC)/analyze.c 
	$(CC) $(CFLAGS) -c  $(SRC)/analyze.c -o $@

$(OBJ)/pattern.o :  $(SRC)/pattern.c 
	$(CC) $(CFLAGS) -c  $(SRC)/pattern.c -o $@

$(OBJ)/xml2ps.o :  $(SRC)/xml2ps.c 
	$(CC) $(CFLAGS) -c  $(SRC)/xml2ps.c -o $@

$(OBJ)/multiple.o :  $(SRC)/multiple.c
	$(CC) $(CFLAGS) -c  $(SRC)/multiple.c -o $@

$(OBJ)/statistics.o :  $(SRC)/statistics.c 
	$(CC) $(CFLAGS) -c  $(SRC)/statistics.c -o $@


clean:
	@rm -f $(OBJ)/*.o
	@rm -f $(ALLTARGETS)

export:
	mkdir -p $(EXPORT_DIR)
	@cd $(EXPORT_DIR); mkdir -p $(INCL)
	$(EXPORT) $(EXPORT_LIST) $(HFILES) $(EXPORT_DIR)/$(INCL)
	@cd $(EXPORT_DIR); mkdir -p $(SRC)
	$(EXPORT) $(EXPORT_LIST) $(SRCFILES) $(TARGETSRC) $(EXPORT_DIR)/$(SRC)
	@cd $(EXPORT_DIR); mkdir -p $(BIN)
	@cd $(EXPORT_DIR); mkdir -p $(LIB)
	@cd $(EXPORT_DIR); mkdir -p $(OBJ)
	@cd $(EXPORT_DIR); mkdir -p $(TEST)
	@cp $(TESTFILES) $(EXPORT_DIR)/$(TEST)
	@cp Makefile $(EXPORT_DIR)




