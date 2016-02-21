SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDLIBS = -lm

all : hinv2

install : hinv2 hinv2.o hinv2.h 
	ar rv $(libdir)/libds.a hinv2.o
	cp hinv2.h $(includedir)
	cp hinv2 $(bindir)

hinv2 : hinv2_main.o hinv2.o $(libdir)/libds.a

hinv2.o : $(includedir)/rtbp.h

clean : 
	rm hinv2 hinv2_main.o hinv2.o
