SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : porbits portbp portbpsym

install : porbits portbp portbpsym portbpsym.o portbpsym.h
	cp porbits portbp portbpsym $(bindir)
	ar rv $(libdir)/libds.a portbp.o portbpsym.o
	cp portbp.h portbpsym.h $(includedir)

porbits : porbits.o portbpsym.o $(libdir)/libds.a

porbits.o : $(includedir)/initcond.h $(includedir)/prtbp_2d.h portbpsym.h 

portbp : portbp_main.o portbp.o $(libdir)/libds.a

portbp_main.o : portbp.h

portbp.o : $(includedir)/prtbp_2d.h $(includedir)/rtbp.h

portbpsym : portbpsym_main.o portbpsym.o $(libdir)/libds.a

portbpsym_main.o : portbpsym.h

portbpsym.o : $(includedir)/prtbp_2d.h $(includedir)/rtbp.h

clean : 
	rm porbits porbits.o \
	   portbp portbp_main.o portbp.o \
	   portbpsym portbpsym_main.o portbpsym.o
