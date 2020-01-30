SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3 -L$(libdir)
LDLIBS = -lgsl -lgslcblas -lds -lm

RESULTS = $(UNSTRESULTS) $(STRESULTS)

all : prtbp_nl_2d prtbp_nl.o prtbp_nl_2d_module.o

install : prtbp_nl_2d prtbp_nl.o prtbp_nl_2d_module.o prtbp_nl.h \
   prtbp_nl_2d_module.h
	ar rv $(libdir)/libds.a prtbp_nl.o prtbp_nl_2d_module.o
	cp prtbp_nl.h prtbp_nl_2d_module.h $(includedir)
	cp prtbp_nl_2d $(bindir)

prtbp_nl.o : $(includedir)/frtbp.h $(includedir)/rtbp.h

prtbp_nl_2d : prtbp_nl_2d.o prtbp_nl_2d_module.o prtbp_nl.o

prtbp_nl_2d.o : prtbp_nl_2d_module.h

prtbp_nl_2d_module.o : $(includedir)/rtbp.h $(includedir)/hinv.h prtbp_nl.h 

clean : 
	rm prtbp_nl_2d prtbp_nl_2d.o prtbp_nl.o \
	   prtbp_nl_2d_module.o
