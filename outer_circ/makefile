SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lds -lgsl -lgslcblas

all : outer_circ outer_circ_test outer_circ_test_neg

install : outer_circ outer_circ_test outer_circ_test_neg
	#ar rv $(libdir)/libds.a inner_ell.o
	#cp inner_ell.h $(includedir)
	cp outer_circ outer_circ_test outer_circ_test_neg $(bindir)

outer_circ : outer_circ.o $(libdir)/libds.a
#	$(CC) -o inner_circ $(LDLIBS) $(CFLAGS) inner_circ_main.o inner_circ.o rtbp.o

outer_circ.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

outer_circ_test : outer_circ_test.o $(libdir)/libds.a
#	$(CC) -o inner_circ $(LDLIBS) $(CFLAGS) inner_circ_main.o inner_circ.o rtbp.o

outer_circ_test.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

outer_circ_test_neg : outer_circ_test_neg.o $(libdir)/libds.a
#	$(CC) -o inner_circ $(LDLIBS) $(CFLAGS) inner_circ_main.o inner_circ.o rtbp.o

outer_circ_test_neg.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

clean : 
	rm outer_circ outer_circ_test  outer_circ_test_neg \
	   outer_circ.o outer_circ_test.o outer_circ_test_neg.o
