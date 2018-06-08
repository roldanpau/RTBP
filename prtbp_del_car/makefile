SHELL = /bin/sh
prefix = /usr/local
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

PROGS = prtbp_del_car prtbp_del_car_inv
OBJS = prtbp_del_car.o
HEADERS = prtbp_del_car.h

all : $(PROGS)

install : $(PROGS) $(OBJS) $(HEADERS)
	ar rv $(libdir)/libds.a $(OBJS)
	cp $(HEADERS) $(includedir)
	cp $(PROGS) $(bindir)

results : prtbp_2d prtbp_2d_inv $(RESULTS)

prtbp_del_car : prtbp_del_car_main.o prtbp_del_car.o

prtbp_del_car_inv : prtbp_del_car_inv_main.o prtbp_del_car.o
	$(CC) prtbp_del_car_inv_main.o prtbp_del_car.o $(LOADLIBES) $(LDLIBS) -o $@

prtbp_del_car_main.o : $(includedir)/rtbp.h prtbp_del_car.h

prtbp_del_car_inv_main.o : $(includedir)/rtbp.h prtbp_del_car.h

prtbp_del_car.o : $(includedir)/frtbp.h $(includedir)/rtbp.h

unstmfld_it0.res: unstmfld_it0.dat prtbp_2d
	./prtbp_2d < $< > $@

stmfld_it0_neg.res: stmfld_it0_neg.dat prtbp_2d_inv
	./prtbp_2d_inv < $< > $@

clean : 
	rm $(PROGS) $(OBJS) prtbp_del_car_main.o prtbp_del_car_inv_main.o
