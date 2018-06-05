SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDLIBS = -lds -lm -lgsl -lgslcblas

PROGS = inner_ell_stoch

all : $(PROGS)

install : $(PROGS)
	ar rv $(libdir)/libds.a inner_ell_stoch.o
	cp inner_ell_stoch.h $(includedir)
	cp $(PROGS) $(bindir)

inner_ell_stoch : inner_ell_stoch_main.o inner_ell_stoch.o
#	$(CC) -o inner_ell_stoch $(LDLIBS) $(CFLAGS) inner_ell_stoch_main.o inner_ell_stoch.o \
#	$(libdir)/libds.a

inner_ell_stoch_main.o : inner_ell_stoch.h

inner_ell_stoch.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

clean : 
	rm $(PROGS) \
	inner_ell_stoch_main.o inner_ell_stoch.o
