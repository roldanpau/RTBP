SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lgsl -lgslcblas -lds -lm

RESULTS = $(UNSTRESULTS) $(STRESULTS)

UNSTRESULTS = unstmfld_it0.res

STRESULTS = stmfld_it0_neg.res

all : prtbp prtbp_inv prtbp_2d prtbp_2d_inv prtbps_2d

install : prtbp prtbp_inv prtbp_2d prtbp_2d_inv prtbp.o prtbp_2d.o prtbp.h \
   prtbp_2d.h
	ar rv $(libdir)/libds.a prtbp.o prtbp_2d.o
	cp prtbp.h prtbp_2d.h $(includedir)
	cp prtbp prtbp_inv prtbp_2d prtbp_2d_inv $(bindir)

results : prtbp_2d prtbp_2d_inv $(RESULTS)

prtbp : prtbp_main.o prtbp.o

prtbp_main.o : $(includedir)/rtbp.h prtbp.h

prtbp.o : $(includedir)/frtbp.h $(includedir)/rtbp.h

prtbp_inv : prtbp_inv.o prtbp.o 

prtbp_inv.o : $(includedir)/rtbp.h prtbp.h

prtbp_2d : prtbp_2d_main.o prtbp_2d.o prtbp.o
#	$(CC) -o prtbp_2d $(LDLIBS) $(CFLAGS) prtbp_2d_main.o prtbp_2d.o \
#	prtbp.o

prtbp_2d_main.o : prtbp_2d.h

prtbp_2d.o : $(includedir)/rtbp.h $(includedir)/hinv.h prtbp.h 

prtbp_2d_inv : prtbp_2d_inv.o prtbp_2d.o prtbp.o

prtbp_2d_inv.o : prtbp_2d.h

prtbps_2d : prtbps_2d.o prtbp_2d.o prtbp.o

prtbps_2d.o : prtbp_2d.h

unstmfld_it0.res: unstmfld_it0.dat prtbp_2d
	./prtbp_2d < $< > $@

stmfld_it0_neg.res: stmfld_it0_neg.dat prtbp_2d_inv
	./prtbp_2d_inv < $< > $@

clean : 
	rm prtbp prtbp_main.o prtbp.o \
	   prtbp_2d prtbp_2d_main.o prtbp_2d.o \
	   prtbps_2d prtbps_2d.o \
	   prtbp_inv prtbp_inv.o \
	   prtbp_2d_inv prtbp_2d_inv.o
