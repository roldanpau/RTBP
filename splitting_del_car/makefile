SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : niterates splitting_del_car

install : splitting_del_car
	cp splitting_del_car $(bindir)
	ar rv $(libdir)/libds.a splitting_del_car.o
	cp splitting_del_car.h $(includedir)

niterates : niterates.o $(libdir)/libds.a

splitting_del_car : splitting_del_car_main.o splitting_del_car.o $(libdir)/libds.a

splitting_del_car.o : $(includedir)/prtbp_2d.h $(includedir)/dprtbp_2d.h

clean : 
	rm niterates splitting_del_car niterates.o splitting_del_car_main.o splitting_del_car.o
