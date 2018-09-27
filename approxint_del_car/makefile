SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

all : approxint_del_car

install : approxint_del_car
	ar rv $(libdir)/libds.a approxint_del_car.o
	cp approxint_del_car.h $(includedir)
	cp approxint_del_car $(bindir)

approxint_del_car : approxint_del_car_main.o approxint_del_car.o

approxint_del_car.o : $(includedir)/prtbp_2d.h

clean : 
	rm approxint_del_car approxint_del_car_main.o approxint_del_car.o
