SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDFLAGS = -O3 #-g
LDLIBS = -lds -lgsl -lgslcblas -lm

all : hinvdel

install : hinvdel hinvdel.o hinvdel.h 
	ar rv $(libdir)/libds.a hinvdel.o
	cp hinvdel.h $(includedir)
	cp hinvdel $(bindir)

hinvdel : hinvdel_main.o hinvdel.o

hinvdel.o : $(includedir)/rtbpdel.h

clean : 
	rm hinvdel hinvdel_main.o hinvdel.o
