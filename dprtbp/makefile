SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

all : dprtbp dprtbp_test dprtbp_2d dprtbp_2d_test

install : dprtbp.o dprtbp.h dprtbp dprtbp_2d.o dprtbp_2d.h dprtbp_2d
	ar rv $(libdir)/libds.a dprtbp.o dprtbp_2d.o
	cp dprtbp.h dprtbp_2d.h $(includedir)
	cp dprtbp dprtbp_2d $(bindir)

dprtbp : dprtbp_main.o dprtbp.o

dprtbp_test : dprtbp_test.o dprtbp.o

dprtbp_test.o : $(includedir)/rtbp.h $(includedir)/prtbp.h dprtbp.h

dprtbp_main.o : $(includedir)/rtbp.h

dprtbp.o : $(includedir)/prtbp.h $(includedir)/frtbp.h $(includedir)/rtbp.h

dprtbp_2d : dprtbp_2d_main.o dprtbp_2d.o dprtbp.o

dprtbp_2d_main.o : dprtbp_2d.h

dprtbp_2d.o : $(includedir)/rtbp.h $(includedir)/prtbp.h \
   $(includedir)/hinv.h dprtbp.h dprtbp_2d.h

dprtbp_2d_test : dprtbp_2d_test.o dprtbp_2d.o dprtbp.o

dprtbp_2d_test.o : $(includedir)/prtbp_2d.h dprtbp_2d.h

clean : 
	rm dprtbp dprtbp_main.o dprtbp.o \
	   dprtbp_2d dprtbp_2d_main.o dprtbp_2d.o \
	   dprtbp_test dprtbp_test.o \
	   dprtbp_2d_test dprtbp_2d_test.o
