SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lds -lgsl -lgslcblas

all : outer_ell B_f B_b

install : outer_ell B_f B_b
#	ar rv $(libdir)/libds.a inner_circ.o outer_ell.o
#	cp inner_circ.h outer_ell.h $(includedir)
	cp outer_ell B_f B_b $(bindir)

re_integrand_B : re_integrand_B.o outer_ell.o $(libdir)/libds.a
	$(CC) -o re_integrand_B $(LDLIBS) $(CFLAGS) re_integrand_B.o outer_ell.o \
	$(libdir)/libds.a

outer_ell : outer_ell_main.o outer_ell.o $(libdir)/libds.a
	$(CC) -o outer_ell $(LDLIBS) $(CFLAGS) outer_ell_main.o outer_ell.o \
	$(libdir)/libds.a

B_f : B_f.o

B_b : B_b.o

outer_ell_main.o : $(includedir)/rtbp.h outer_ell.h

outer_ell.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h \
	$(includedir)/inner_ell.h 

clean : 
	rm re_integrand_B outer_ell B_f B_b \
	   re_integrand_B.o outer_ell_main.o outer_ell.o B_f.o B_b.o
