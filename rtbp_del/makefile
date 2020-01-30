SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDFLAGS = -O3 #-g
LDLIBS = -lm -lgsl -lgslcblas

all : rtbpdel

install : rtbpdel rtbpdel.o rtbpdel.h
	ar rv $(libdir)/libds.a rtbpdel.o
	cp rtbpdel.h $(includedir)
	cp rtbpdel $(bindir)

rtbpdel : rtbpdel_main.o rtbpdel.o

rtbpdel_main.o : rtbpdel.h

rtbpdel.o : rtbpdel.h

clean : 
	rm rtbpdel rtbpdel_main.o rtbpdel.o
