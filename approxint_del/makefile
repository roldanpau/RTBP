SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : approxintdel

install : approxintdel
	cp approxintdel $(bindir)
	ar rv $(libdir)/libds.a approxintdel.o
	cp approxintdel.h $(includedir)

approxintdel : approxintdel_main.o approxintdel.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

approxintdel_main.o : 

approxintdel.o : $(includedir)/prtbpdel_2d.h

clean : 
	rm approxintdel approxintdel_main.o approxintdel.o
