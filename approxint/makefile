SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : approxint

install : approxint
	cp approxint $(bindir)
	ar rv $(libdir)/libds.a approxint.o
	cp approxint.h $(includedir)

approxint : approxint_main.o approxint.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

approxint_main.o : 

approxint.o : $(includedir)/prtbp_2d.h

clean : 
	rm approxint approxint_main.o approxint.o
