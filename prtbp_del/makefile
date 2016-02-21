SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

RESULTS = $(UNSTRESULTS) $(STRESULTS)

UNSTRESULTS = 

STRESULTS = 

PROGS = prtbpdel prtbpdel_inv \
	prtbpdel_2d prtbpdel_2d_inv \
	sec1sec2del sec1sec2del_inv

all : $(PROGS)

install : $(PROGS) \
   prtbpdel.o prtbpdel.h \
   prtbpdel_2d.o prtbpdel_2d.h
	ar rv $(libdir)/libds.a prtbpdel.o prtbpdel_2d.o
	cp prtbpdel.h prtbpdel_2d.h $(includedir)
	cp $(PROGS) $(bindir)

results : $(PROGS) $(RESULTS)
	
prtbpdel : prtbpdel_main.o prtbpdel.o $(libdir)/libds.a

prtbpdel_main.o : $(includedir)/rtbp.h prtbpdel.h

prtbpdel_inv : prtbpdel_inv_main.o prtbpdel.o $(libdir)/libds.a
	$(CC) -o prtbpdel_inv $(LDLIBS) $(CFLAGS) prtbpdel_inv_main.o \
	   prtbpdel.o $(libdir)/libds.a

prtbpdel_inv_main.o : $(includedir)/rtbp.h prtbpdel.h

prtbpdel.o : $(includedir)/frtbpdel.h $(includedir)/rtbp.h

prtbpdel_2d : prtbpdel_2d_main.o prtbpdel_2d.o prtbpdel.o $(libdir)/libds.a

sec1sec2del: sec1sec2del_main.o prtbpdel_2d.o prtbpdel.o $(libdir)/libds.a
	$(CC) -o $@ $(LDLIBS) $(CFLAGS) sec1sec2del_main.o prtbpdel_2d.o \
	prtbpdel.o $(libdir)/libds.a

sec1sec2del_inv: sec1sec2del_inv_main.o prtbpdel_2d.o prtbpdel.o $(libdir)/libds.a
	$(CC) -o $@ $(LDLIBS) $(CFLAGS) sec1sec2del_inv_main.o prtbpdel_2d.o \
	prtbpdel.o $(libdir)/libds.a

prtbpdel_2d_main.o : prtbpdel_2d.h

sec1sec2del_main.o : prtbpdel.h prtbpdel_2d.h

sec1sec2del_inv_main.o : prtbpdel.h prtbpdel_2d.h

prtbpdel_2d_inv : prtbpdel_2d_inv_main.o prtbpdel_2d.o prtbpdel.o $(libdir)/libds.a
	$(CC) -o prtbpdel_2d_inv $(LDLIBS) $(CFLAGS) prtbpdel_2d_inv_main.o \
	   prtbpdel_2d.o prtbpdel.o $(libdir)/libds.a

prtbpdel_2d_inv_main.o : prtbpdel_2d.h

prtbpdel_2d.o : $(includedir)/hinvdel.h $(includedir)/rtbp.h prtbpdel.h

unstmfld_it0.res: unstmfld_it0.dat prtbpdel
	./prtbpdel < $< > $@

stmfld_it0_neg.res: stmfld_it0_neg.dat prtbpdel_inv
	./prtbpdel_inv < $< > $@

clean : 
	rm $(PROGS) \
	   prtbpdel_main.o prtbpdel.o \
	   prtbpdel_inv_main.o \
	   prtbpdel_2d_main.o prtbpdel_2d.o \
	   prtbpdel_2d_inv_main.o \
           sec1sec2del_main.o \
           sec1sec2del_inv_main.o \
	   $(RESULTS)
