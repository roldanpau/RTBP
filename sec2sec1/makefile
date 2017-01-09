SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
#LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

PROGS = sec2sec1

all : $(PROGS)

install : $(PROGS) sec2sec1_module.o sec2sec1_module.h
	cp $(PROGS) $(bindir)
	ar rv $(libdir)/libds.a sec2sec1_module.o
	cp sec2sec1_module.h $(includedir)

sec2sec1: sec2sec1_module.o $(libdir)/libds.a

sec2sec1_inv: sec2sec1_module.o $(libdir)/libds.a

sec2sec1_module.o : $(includedir)/rtbp.h $(includedir)/hinv.h

clean : 
	rm $(PROGS) \
	   sec2sec1_module.o \
