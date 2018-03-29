SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lds -lm

all : frtbp

install : frtbp frtbp.o frtbp.h
	ar rv $(libdir)/libds.a frtbp.o
	cp frtbp.h $(includedir)
	cp frtbp $(bindir)

frtbp : frtbp.o frtbp_main.o
#	$(CC) -o frtbp $(LDLIBS) $(CFLAGS) frtbp_main.o frtbp.o rtbp.o

frtbp_main.o : frtbp.h

frtbp.o : $(includedir)/rtbp.h

clean : 
	rm frtbp frtbp_main.o frtbp.o
