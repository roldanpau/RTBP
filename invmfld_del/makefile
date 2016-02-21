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
	  stmfld.res stmfld_neg.res \
	  unstmfld_it0.res unstmfld_it0_neg.res \
	  stmfld_it0.res stmfld_it0_neg.res \
	  unstmfld_it2.res stmfld_it3.res

all : invmflddel

install : invmflddel
	cp invmflddel $(bindir)

results: invmflddel $(RESULTS)

invmflddel : invmflddel.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

invmflddel.o : $(includedir)/disc.h

%.res: %.dat invmflddel
	./invmflddel < $< > $@

clean : 
	rm invmflddel $(RESULTS) invmflddel.o
