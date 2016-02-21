SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDLIBS = -lm -lds -lgsl -lgslcblas

all : re_integrand_A inner_ell

install : inner_ell
	ar rv $(libdir)/libds.a inner_ell.o
	cp inner_ell.h $(includedir)
	cp inner_ell $(bindir)

re_integrand_A : re_integrand_A.o $(libdir)/libds.a
	$(CC) -o re_integrand_A $(LDLIBS) $(CFLAGS) re_integrand_A.o \
	$(libdir)/libds.a

re_integrand_A.o : inner_ell.h

inner_ell : inner_ell_main.o inner_ell.o $(libdir)/libds.a
	$(CC) -o inner_ell $(LDLIBS) $(CFLAGS) inner_ell_main.o inner_ell.o \
	$(libdir)/libds.a

inner_ell_main.o : inner_ell.h

inner_ell.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

clean : 
	rm re_integrand_A inner_ell \
	re_integrand_A.o inner_ell_main.o inner_ell.o
