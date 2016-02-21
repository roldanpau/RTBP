SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

PROGS = portbpdel porbitsdel

all : $(PROGS)

install : portbpdel portbpdel.o portbpdel.h
	cp portbpdel $(bindir)
	ar rv $(libdir)/libds.a portbpdel.o 
	cp portbpdel.h $(includedir)

porbitsdel : porbitsdel.o portbpdel.o $(libdir)/libds.a

porbitsdel.o : $(includedir)/rtbp.h $(includedir)/prtbpdel_2d.h portbpdel.h 

portbpdel : portbpdel_main.o portbpdel.o $(libdir)/libds.a

portbpdel_main.o : portbpdel.h

portbpdel.o : $(includedir)/prtbpdel_2d.h $(includedir)/dprtbpdel_2d.h \
   $(includedir)/rtbp.h

clean : 
	rm $(PROGS) portbpdel_main.o portbpdel.o porbitsdel.o
