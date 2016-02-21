SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm

all : hinv

install : hinv hinv.o hinv.h 
	ar rv $(libdir)/libds.a hinv.o
	cp hinv.h $(includedir)
	cp hinv $(bindir)

hinv : hinv_main.o hinv.o

hinv.o : $(includedir)/rtbp.h

clean : 
	rm hinv hinv_main.o hinv.o
