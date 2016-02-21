SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib

# WARNING! Do not activate optimization flag -O3, because it makes cardel2_2d
# numerically unstable, for some reason!
# NOTE: Now I'm not so sure it's due to optimization, so I let it be...
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

RESULTS = \
	  unstmfld_it0.res unstmfld_it0_neg.res \
	  stmfld_it0.res stmfld_it0_neg.res \
	  unstmfld_it2.res unstmfld_it2_neg.res \
	  stmfld_it2.res stmfld_it2_neg.res \
          unstmfld_it3_neg.res stmfld_it3_neg.res


all : cardel2s_2d cardel2_2d

install : cardel2s_2d cardel2_2d cardel2_2d.o cardel2_2d.h
	cp cardel2s_2d cardel2_2d $(bindir)
	ar rv $(libdir)/libds.a cardel2_2d.o
	cp cardel2_2d.h $(includedir)

results : cardel2_2d $(RESULTS)

cardel2s_2d : cardel2s_2d_main.o cardel2_2d.o $(libdir)/libds.a
	$(CC) -o cardel2s_2d $(LDLIBS) $(CFLAGS) cardel2s_2d_main.o \
	   cardel2_2d.o $(libdir)/libds.a

cardel2_2d : cardel2_2d_main.o cardel2_2d.o $(libdir)/libds.a
	$(CC) -o cardel2_2d $(LDLIBS) $(CFLAGS) cardel2_2d_main.o \
	   cardel2_2d.o $(libdir)/libds.a

cardel2s_2d_main.o : $(includedir)/rtbp.h cardel2_2d.h

cardel2_2d_main.o : $(includedir)/rtbp.h cardel2_2d.h

cardel2_2d.o : $(includedir)/rtbp.h $(includedir)/hinv2.h $(includedir)/cardel.h

%.res: %.dat cardel2_2d
	./cardel2_2d < $< > $@

clean : 
	rm cardel2s_2d cardel2_2d \
	   cardel2s_2d_main.o cardel2_2d_main.o \
	   cardel2_2d.o \
	   $(RESULTS)
