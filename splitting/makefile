SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : splitting

install : splitting
	cp splitting $(bindir)
	ar rv $(libdir)/libds.a splitting.o
	cp splitting.h $(includedir)

splitting : splitting_main.o splitting.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

splitting_main.o : 

splitting.o : $(includedir)/prtbp_2d.h $(includedir)/dprtbp_2d.h

clean : 
	rm splitting splitting_main.o splitting.o
