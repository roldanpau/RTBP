SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

all : frtbpred 

install : frtbpred frtbpred.o rtbpred.o frtbpred.h rtbpred.h
	ar rv $(libdir)/libds.a frtbpred.o rtbpred.o
	cp frtbpred.h rtbpred.h $(includedir)
	cp frtbpred $(bindir)

frtbpred : frtbpred_main.o frtbpred.o rtbpred.o 
#	$(CC) -o frtbpred $(LDLIBS) $(CFLAGS) frtbpred_main.o frtbpred.o rtbpred.o

frtbpred_main.o : rtbpred.h

frtbpred.o : $(includedir)/rtbp.h

rtbpred.o : rtbpred.h

clean : 
	rm rtbpred.o
	rm frtbpred frtbpred_main.o frtbpred.o
