SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : splittingdel

install : splittingdel
	cp splittingdel $(bindir)
	ar rv $(libdir)/libds.a splittingdel.o
	cp splittingdel.h $(includedir)

splittingdel : splittingdel_main.o splittingdel.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

splittingdel_main.o : 

splittingdel.o : $(includedir)/prtbpdel_2d.h $(includedir)/dprtbpdel_2d.h

clean : 
	rm splittingdel splittingdel_main.o splittingdel.o
