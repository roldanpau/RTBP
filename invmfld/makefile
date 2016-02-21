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

all : invmfld

install : invmfld
	cp invmfld $(bindir)
	ar rv $(libdir)/libds.a disc.o
	cp disc.h $(includedir)

results: invmfld $(RESULTS)

invmfld : invmfld.o disc.o $(libdir)/libds.a

invmfld.o : $(includedir)/prtbp_2d.h

%.res: %.dat invmfld
	./invmfld < $< > $@

clean : 
	rm invmfld $(RESULTS) invmfld.o disc.o
