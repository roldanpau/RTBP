SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : intersecdel

install : intersecdel
	cp intersecdel $(bindir)
	ar rv $(libdir)/libds.a intersecdel.o
	cp intersecdel.h $(includedir)

intersecdel : intersecdel_main.o intersecdel.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

intersecdel_main.o : 

intersecdel.o : $(includedir)/prtbpdel_2d.h

clean : 
	rm intersecdel intersecdel_main.o intersecdel.o
