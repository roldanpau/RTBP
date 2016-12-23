SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

RESULTS = \
	  unstmfld.res unstmfld_neg.res \
	  stmfld.res stmfld_neg.res

all : invmfld_del_car

install : invmfld_del_car
	cp invmfld_del_car $(bindir)
	ar rv $(libdir)/libds.a lift.o
	cp lift.h $(includedir)

results: invmfld_del_car $(RESULTS)

invmfld_del_car : invmfld_del_car.o lift.o $(libdir)/libds.a

invmfld_del_car.o : $(includedir)/prtbp_2d.h

%.res: %.dat invmfld_del_car
	./invmfld_del_car < $< > $@

clean : 
	rm invmfld_del_car $(RESULTS) invmfld_del_car.o lift.o
