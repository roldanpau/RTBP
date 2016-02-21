SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

RESULTS = $(UNSTRESULTS) $(STRESULTS)

UNSTRESULTS = \
	unstmfld_it1.res unstmfld_it1_neg.res \
	unstmfld_it2.res unstmfld_it2_neg.res \
	unstmfld_it3.res unstmfld_it3_neg.res \
	unstmfld_it4.res unstmfld_it4_neg.res \
	unstmfld_it5.res unstmfld_it5_neg.res \
	unstmfld_it6.res unstmfld_it6_neg.res \
	unstmfld_it7.res unstmfld_it7_neg.res

STRESULTS = \
	stmfld_it1.res stmfld_it1_neg.res \
	stmfld_it2.res stmfld_it2_neg.res \
	stmfld_it3.res stmfld_it3_neg.res \
	stmfld_it4.res stmfld_it4_neg.res \
	stmfld_it5.res stmfld_it5_neg.res \
	stmfld_it6.res stmfld_it6_neg.res \
	stmfld_it7.res stmfld_it7_neg.res

all : orbitpdel orbitpdel_bwd

#install : orbitpdel orbitpdel_bwd orbitpdel.o orbitpdel.h
#	cp orbitpdel orbitpdel_bwd $(bindir)
#	ar rv $(libdir)/libds.a orbitpdel.o
#	cp orbitpdel.h $(includedir)

results: orbitpdel orbitpdel_bwd $(RESULTS)

orbitpdel : orbitpdel_main.o orbitpdel.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

orbitpdel_bwd : orbitpdel_bwd_main.o orbitpdel.o $(libdir)/libds.a
	$(CC) -o orbitpdel_bwd $(LDLIBS) $(CFLAGS) orbitpdel_bwd_main.o \
	   orbitpdel.o $(libdir)/libds.a

orbitpdel_main.o : $(includedir)/prtbpdel.h orbitpdel.h

orbitpdel_bwd_main.o : $(includedir)/prtbpdel.h orbitpdel.h

orbitpdel.o : 

unstmfld_it1.res: unstmfld_it1.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it1_neg.res: unstmfld_it1_neg.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it2.res: unstmfld_it2.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it2_neg.res: unstmfld_it2_neg.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it3.res: unstmfld_it3.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it3_neg.res: unstmfld_it3_neg.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it4.res: unstmfld_it4.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it4_neg.res: unstmfld_it4_neg.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it5.res: unstmfld_it5.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it5_neg.res: unstmfld_it5_neg.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it6.res: unstmfld_it6.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it6_neg.res: unstmfld_it6_neg.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it7.res: unstmfld_it7.dat orbitpdel
	./orbitpdel < $< > $@

unstmfld_it7_neg.res: unstmfld_it7_neg.dat orbitpdel
	./orbitpdel < $< > $@

stmfld_it1.res: stmfld_it1.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it1_neg.res: stmfld_it1_neg.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it2.res: stmfld_it2.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it2_neg.res: stmfld_it2_neg.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it3.res: stmfld_it3.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it3_neg.res: stmfld_it3_neg.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it4.res: stmfld_it4.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it4_neg.res: stmfld_it4_neg.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it5.res: stmfld_it5.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it5_neg.res: stmfld_it5_neg.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it6.res: stmfld_it6.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it6_neg.res: stmfld_it6_neg.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it7.res: stmfld_it7.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

stmfld_it7_neg.res: stmfld_it7_neg.dat orbitpdel_bwd
	./orbitpdel_bwd < $< > $@

clean : 
	rm orbitpdel orbitpdel_bwd \
	   $(RESULTS) \
	   orbitpdel_main.o orbitpdel_bwd_main.o orbitpdel.o
