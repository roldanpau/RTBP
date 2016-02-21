SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : Ldeviation

install : Ldeviation Ldeviation.o Ldeviation.h
	ar rv $(libdir)/libds.a Ldeviation.o 
	cp Ldeviation.h $(includedir)
	cp Ldeviation $(bindir)

Ldeviation : Ldeviation_main.o Ldeviation.o $(libdir)/libds.a

Ldeviation_main.o : $(includedir)/rtbp.h Ldeviation.h

Ldeviation.o : $(includedir)/frtbpdel.h $(includedir)/rtbp.h

clean : 
	rm Ldeviation Ldeviation_main.o Ldeviation.o
