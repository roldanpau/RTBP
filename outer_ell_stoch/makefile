SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDLIBS = -lds -lm -lgsl -lgslcblas

PROGS = outer_ell_stoch B_j

all : $(PROGS)

install : $(PROGS)
	#ar rv $(libdir)/libds.a outer_ell_stoch.o
	#cp outer_ell_stoch.h $(includedir)
	cp $(PROGS) $(bindir)

outer_ell_stoch : outer_ell_stoch_main.o outer_ell_stoch.o
#	$(CC) -o outer_ell_stoch $(LDLIBS) $(CFLAGS) outer_ell_stoch_main.o outer_ell_stoch.o \
#	$(libdir)/libds.a

B_j : B_j.o

outer_ell_stoch_main.o : outer_ell_stoch.h

outer_ell_stoch.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

clean : 
	rm $(PROGS) \
	outer_ell_stoch_main.o outer_ell_stoch.o B_j.o
