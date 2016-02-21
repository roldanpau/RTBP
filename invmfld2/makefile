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

all : invmfld2

install : invmfld2
	cp invmfld2 $(bindir)
	#ar rv $(libdir)/libds.a disc.o
	#cp disc.h $(includedir)

results: invmfld2 $(RESULTS)

invmfld2 : invmfld2.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

invmfld2.o : $(includedir)/orbitp.h

%.res: %.dat invmfld2
	./invmfld2 < $< > $@

clean : 
	rm invmfld2 $(RESULTS) invmfld2.o
