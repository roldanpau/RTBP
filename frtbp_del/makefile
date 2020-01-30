SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDFLAGS = -O3 #-g
LDLIBS = -lgsl -lgslcblas -lds -lm

all : frtbpdel 

install : frtbpdel frtbpdel.o frtbpdel.h 
	ar rv $(libdir)/libds.a frtbpdel.o 
	cp frtbpdel.h $(includedir)
	cp frtbpdel $(bindir)

frtbpdel : frtbpdel_main.o frtbpdel.o
#	$(CC) -o frtbpdel $(LDLIBS) $(CFLAGS) frtbpdel_main.o frtbpdel.o rtbpdel.o

frtbpdel_main.o : $(includedir)/rtbpdel.h

frtbpdel.o : $(includedir)/rtbp.h

clean : 
	rm frtbpdel frtbpdel_main.o frtbpdel.o
