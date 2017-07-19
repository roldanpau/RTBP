SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm

all : lift.o lift.h hinv

install : hinv hinv.o hinv.h lift.o lift.h
	ar rv $(libdir)/libds.a hinv.o lift.o
	cp hinv.h lift.h $(includedir)
	cp hinv $(bindir)

hinv : hinv_main.o hinv.o

hinv.o : $(includedir)/rtbp.h

lift.o : $(includedir)/rtbp.h $(includedir)/section.h

clean : 
	rm hinv hinv_main.o hinv.o lift.o
