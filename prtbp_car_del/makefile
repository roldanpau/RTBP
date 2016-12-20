SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

RESULTS = $(UNSTRESULTS) $(STRESULTS)

UNSTRESULTS = unstmfld_it0.res

STRESULTS = stmfld_it0_neg.res

all : prtbp_car_del prtbp_car_del_inv 

install : prtbp_car_del prtbp_car_del_inv prtbp_car_del.o prtbp_car_del.h
	ar rv $(libdir)/libds.a prtbp_car_del.o
	cp prtbp_car_del.h $(includedir)
	cp prtbp_car_del prtbp_car_del_inv $(bindir)

results : prtbp_2d prtbp_2d_inv $(RESULTS)

prtbp_car_del : prtbp_car_del_main.o prtbp_car_del.o $(libdir)/libds.a

prtbp_car_del_main.o : $(includedir)/rtbp.h prtbp_car_del.h

prtbp_car_del.o : $(includedir)/frtbp.h $(includedir)/rtbp.h

prtbp_car_del_inv : prtbp_car_del_inv.o prtbp_car_del.o $(libdir)/libds.a

prtbp_car_del_inv.o : $(includedir)/rtbp.h prtbp_car_del.h

prtbp_2d : prtbp_2d_main.o prtbp_2d.o prtbp.o $(libdir)/libds.a
#	$(CC) -o prtbp_2d $(LDLIBS) $(CFLAGS) prtbp_2d_main.o prtbp_2d.o \
#	prtbp.o

prtbp_2d_main.o : prtbp_2d.h

prtbp_2d.o : $(includedir)/rtbp.h $(includedir)/hinv.h prtbp.h 

prtbp_2d_inv : prtbp_2d_inv.o prtbp_2d.o prtbp.o $(libdir)/libds.a

prtbp_2d_inv.o : prtbp_2d.h

prtbps_2d : prtbps_2d.o prtbp_2d.o prtbp.o $(libdir)/libds.a

prtbps_2d.o : prtbp_2d.h

unstmfld_it0.res: unstmfld_it0.dat prtbp_2d
	./prtbp_2d < $< > $@

stmfld_it0_neg.res: stmfld_it0_neg.dat prtbp_2d_inv
	./prtbp_2d_inv < $< > $@

clean : 
	rm prtbp_car_del prtbp_car_del_main.o prtbp_car_del.o \
	   prtbp_2d prtbp_2d_main.o prtbp_2d.o \
	   prtbps_2d prtbps_2d.o \
	   prtbp_car_del_inv prtbp_car_del_inv.o \
	   prtbp_2d_inv prtbp_2d_inv.o \
	   $(RESULTS)
