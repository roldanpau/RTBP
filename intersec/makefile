SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lgsl -lgslcblas -lds -lm

all : intersec

install : intersec
	cp intersec $(bindir)
	ar rv $(libdir)/libds.a intersec.o
	cp intersec.h $(includedir)

intersec : intersec_main.o intersec.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

intersec_main.o : 

intersec.o : $(includedir)/prtbp_2d.h

clean : 
	rm intersec intersec_main.o intersec.o
