SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

all : intersec_del_car

install : intersec_del_car
	cp intersec_del_car $(bindir)
	ar rv $(libdir)/libds.a intersec_del_car.o
	cp intersec_del_car.h $(includedir)

intersec_del_car : intersec_del_car_main.o intersec_del_car.o
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

intersec_del_car.o : $(includedir)/prtbp_2d.h

clean : 
	rm intersec_del_car intersec_del_car_main.o intersec_del_car.o
