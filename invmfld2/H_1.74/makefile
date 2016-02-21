SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

RESULTS = \
	  unstmfld.res unstmfld_neg.res \
	  stmfld.res stmfld_neg.res

results: $(RESULTS)

%.res: %.dat
	invmfld2 < $< > $@

clean : 
	rm $(RESULTS)
