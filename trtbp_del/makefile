SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
LDLIBS = -lm -lgsl -lgslcblas -lds

all : trtbpdel

install : trtbpdel
	cp trtbpdel $(bindir)

trtbpdel : $(libdir)/libds.a
#	cc -o trtbpdel -g -lm -lgsl -lgslcblas trtbpdel.o rtbp.o

#frtbp : frtbp_main.o frtbp.o libds.a(rtbp.o)
#	$(CC) -o frtbp $(LDLIBS) $(CFLAGS) frtbp_main.o frtbp.o rtbp.o

trtbpdel.o : $(includedir)/frtbpdel.h

clean : 
	rm trtbpdel trtbpdel.o
