SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDFLAGS = -O3 #-g
LDLIBS = -lm -lds -lgsl -lgslcblas

all : inner_circ

install : inner_circ
	ar rv $(libdir)/libds.a inner_circ.o
	cp inner_circ.h $(includedir)
	cp inner_circ $(bindir)

inner_circ : inner_circ_main.o inner_circ.o $(libdir)/libds.a
#	$(CC) -o inner_circ $(LDLIBS) $(CFLAGS) inner_circ_main.o inner_circ.o rtbp.o

inner_circ_main.o : inner_circ.h

inner_circ.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

clean : 
	rm inner_circ \
	   inner_circ_main.o inner_circ.o
