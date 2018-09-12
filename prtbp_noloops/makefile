SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -g
LDFLAGS = -g
LDLIBS = -lgsl -lgslcblas -lds -lm

RESULTS = $(UNSTRESULTS) $(STRESULTS)

all : prtbp_nl.o prtbp_nl_2d.o

install : prtbp_nl.o prtbp_nl_2d.o prtbp_nl.h \
   prtbp_nl_2d.h
	ar rv $(libdir)/libds.a prtbp_nl.o prtbp_nl_2d.o
	cp prtbp_nl.h prtbp_nl_2d.h $(includedir)

prtbp_nl.o : $(includedir)/frtbp.h $(includedir)/rtbp.h

prtbp_nl_2d.o : $(includedir)/rtbp.h $(includedir)/hinv.h prtbp_nl.h 

clean : 
	rm prtbp_nl.o \
	   prtbp_nl_2d.o
