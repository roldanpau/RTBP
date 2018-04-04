SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -g
LDFLAGS = -g
LDLIBS = -lgsl -lgslcblas -lds -lm

RESULTS = $(UNSTRESULTS) $(STRESULTS)

UNSTRESULTS = unstmfld_it0.res

STRESULTS = stmfld_it0_neg.res

all : prtbp_del_car

install : prtbp_del_car prtbp_del_car.o prtbp_del_car.h
	ar rv $(libdir)/libds.a prtbp_del_car.o
	cp prtbp_del_car.h $(includedir)
	cp prtbp_del_car $(bindir)

results : prtbp_2d prtbp_2d_inv $(RESULTS)

prtbp_del_car : prtbp_del_car_main.o prtbp_del_car.o

prtbp_del_car_main.o : $(includedir)/rtbp.h prtbp_del_car.h

prtbp_del_car.o : $(includedir)/frtbp.h $(includedir)/rtbp.h

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
	rm prtbp_del_car prtbp_del_car_main.o prtbp_del_car.o 
