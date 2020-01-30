SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

RESULTS = \
	  unstmfld.res unstmfld_neg.res \
	  stmfld.res stmfld_neg.res

all : invmfld

install : invmfld
	cp invmfld $(bindir)
	ar rv $(libdir)/libds.a disc.o
	cp disc.h $(includedir)

results: invmfld $(RESULTS)

invmfld : invmfld.o disc.o

invmfld.o : $(includedir)/prtbp_2d.h

%.res: %.dat invmfld
	./invmfld < $< > $@

clean : 
	rm invmfld invmfld.o disc.o
