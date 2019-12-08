SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -g #-O3
LDFLAGS = -g #-O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : porbits_apo

install : porbits_apo
	cp porbits_apo $(bindir)

porbits_apo : porbits_apo.o $(libdir)/libds.a

porbits_apo.o : $(includedir)/initcond_apo.h $(includedir)/portbpsym.h 

clean : 
	rm porbits_apo porbits_apo.o
