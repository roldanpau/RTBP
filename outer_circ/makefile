SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
#LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

# PRG (2018-03-29)
# What happened with outer_circ_test and outer_circ_test_neg??
# The code files outer_circ_test.c is missing from this dir.
# I don't even remember what whas the purpose of these files.

PROGS = outer_circ outer_circ_stoch outer_circ_stoch_test

all : $(PROGS)

install : $(PROGS)
	ar rv $(libdir)/libds.a \
        outer_circ_stoch_module.o
	cp $(PROGS) $(bindir)
	cp outer_circ_stoch_module.h $(includedir)

outer_circ : outer_circ.o
#	$(CC) -o inner_circ $(LDLIBS) $(CFLAGS) inner_circ_main.o inner_circ.o rtbp.o

outer_circ_stoch : outer_circ_stoch_module.o

outer_circ_stoch_test : outer_circ_stoch_module.o

outer_circ.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

outer_circ_stoch_module.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

clean : 
	rm $(PROGS) \
	   outer_circ.o outer_circ_stoch_module.o
