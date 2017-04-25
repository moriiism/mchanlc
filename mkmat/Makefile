
include ../../mitool/Makefile.in
include ../../mitool/mirootlib/Makefile.in

INCLUDE := $(ROOT_INC) $(MIROOTLIB_INC) $(MILIB_INC) $(MT_INC)
LDLIBS  := $(ROOT_ALLLIBS) $(MIROOTLIB_LIB) $(MILIB_LIB) $(LIB)

TARGET = mkmat  mkmat2

all: $(TARGET) 

mkmat: mkmat.o  arg_mkmat.o
	$(CXXN) $@.o arg_$@.o ${LDLIBS} -o $@

mkmat2: mkmat2.o  arg_mkmat2.o
	$(CXXN) $@.o arg_$@.o ${LDLIBS} -o $@

.SUFFIXES:
.SUFFIXES: .o .c .f .F .cc

.c.o:
	$(CCN) -c $<

.F.o : 
	$(FCN) -c $<

.f.o:
	$(FCN) -c $<

.cc.o:
	${CXXN} -c $<

install: all
	if [ ! -d $(INST_BIN_DIR) ]; then mkdir $(INST_BIN_DIR); fi
	install ${TARGET} ${INST_BIN_DIR}

clean:
	-rm -f *.o

cleaner: clean
	-rm -f ${TARGET} *% *~ core
