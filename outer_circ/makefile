SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lds -lgsl -lgslcblas

# PRG (2018-03-29)
# What happened with outer_circ_test and outer_circ_test_neg??
# The code files outer_circ_test.c is missing from this dir.
# I don't even remember what whas the purpose of these files.

#all : outer_circ outer_circ_test outer_circ_test_neg
all : outer_circ outer_circ_stoch

install : outer_circ outer_circ_stoch
	#ar rv $(libdir)/libds.a inner_ell.o
	#cp inner_ell.h $(includedir)
	cp outer_circ outer_circ_stoch $(bindir)

outer_circ : outer_circ.o $(libdir)/libds.a
#	$(CC) -o inner_circ $(LDLIBS) $(CFLAGS) inner_circ_main.o inner_circ.o rtbp.o

outer_circ_stoch : outer_circ_stoch.o $(libdir)/libds.a

outer_circ.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

outer_circ_stoch.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

clean : 
	rm outer_circ outer_circ_stoch \
	   outer_circ.o outer_circ_stoch.o
