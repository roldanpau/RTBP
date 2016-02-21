SHELL = /bin/sh

hostname := $(shell hostname)
ifeq ($(hostname),pau-PORTEGE-R700)
   prefix = /usr/local
   exec_prefix = $(prefix)
   bindir = $(exec_prefix)/bin
   includedir = $(prefix)/include
   libdir = $(exec_prefix)/lib

   # WARNING! Do not activate optimization flag -O3, because it makes cardel_2d
   # numerically unstable, for some reason!
   # NOTE: Now I'm not so sure it's due to optimization, so I let it be...
   CFLAGS = -O3
   LDFLAGS = -O3
   LDLIBS = -lm -lgsl -lgslcblas -lds
else
   prefix = $(HOME)
   exec_prefix = $(prefix)
   bindir = $(exec_prefix)/bin
   includedir = $(prefix)/include
   libdir = $(exec_prefix)/lib

   # WARNING! Do not activate optimization flag -O3, because it makes cardel_2d
   # numerically unstable, for some reason!
   # NOTE: Now I'm not so sure it's due to optimization, so I let it be...
   CFLAGS = -O3 -I$(GSL_INC) -I$(includedir)
   LDFLAGS = -O3 -L$(GSL_LIB) -L$(libdir) -Wl,-rpath,$(GSL_LIB)
   LDLIBS = -lm -lgsl -lgslcblas -lds
endif

RESULTS = \
	  unstmfld_it0.res unstmfld_it0_neg.res \
	  stmfld_it0.res stmfld_it0_neg.res \
	  unstmfld_it2.res unstmfld_it2_neg.res \
	  stmfld_it2.res stmfld_it2_neg.res \
          unstmfld_it3_neg.res stmfld_it3_neg.res


all : cardels_2d cardel_2d cardel

install : cardels_2d cardel_2d cardel cardel.o cardel_2d.o cardel.h cardel_2d.h
	cp cardels_2d cardel_2d cardel $(bindir)
	ar rv $(libdir)/libds.a cardel.o cardel_2d.o
	cp cardel.h cardel_2d.h $(includedir)

results : cardel_2d $(RESULTS)

cardels_2d : cardels_2d_main.o cardel_2d.o cardel.o $(libdir)/libds.a
	$(CC) -o cardels_2d $(LDLIBS) $(CFLAGS) cardels_2d_main.o \
	   cardel_2d.o cardel.o $(libdir)/libds.a

cardel_2d : cardel_2d_main.o cardel_2d.o cardel.o $(libdir)/libds.a
	$(CC) -o cardel_2d $(LDLIBS) $(CFLAGS) cardel_2d_main.o \
	   cardel_2d.o cardel.o $(libdir)/libds.a

cardel : cardel_main.o cardel.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

cardels_2d_main.o : $(includedir)/rtbp.h cardel_2d.h

cardel_2d_main.o : $(includedir)/rtbp.h cardel_2d.h

cardel_main.o : cardel.h

cardel.o : 

cardel_2d.o : $(includedir)/rtbp.h $(includedir)/hinv.h cardel.h

%.res: %.dat cardel_2d
	./cardel_2d < $< > $@

clean : 
	rm cardels_2d cardel_2d cardel \
	   cardels_2d_main.o cardel_2d_main.o cardel_main.o \
	   cardel_2d.o cardel.o \
	   $(RESULTS)
