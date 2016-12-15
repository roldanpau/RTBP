SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : hyper hypers

install : hyper.o hyper.h hyper hypers
	ar rv $(libdir)/libds.a hyper.o 
	cp hyper.h $(includedir)
	cp hyper hypers $(bindir)

hypers : hypers.o hyper.o $(libdir)/libds.a

hypers.o : 

hyper : hyper_main.o hyper.o $(libdir)/libds.a

hyper_main.o : 

hyper.o : $(includedir)/dprtbp_2d.h

clean : 
	rm hypers hypers.o \
	   hyper hyper_main.o hyper.o 
