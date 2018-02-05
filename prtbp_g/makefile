SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : prtbp_g

install : prtbp_g prtbp_g.o prtbp_g.h
	ar rv $(libdir)/libds.a prtbp_g.o
	cp prtbp_g.h $(includedir)
	cp prtbp_g $(bindir)

prtbp_g : prtbp_g_main.o prtbp_g.o $(libdir)/libds.a

prtbp_g_main.o : $(includedir)/rtbp.h prtbp_g.h

prtbp_g.o : $(includedir)/frtbp.h $(includedir)/rtbp.h

clean : 
	rm prtbp_g prtbp_g_main.o prtbp_g.o
