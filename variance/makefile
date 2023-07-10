SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
#LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

PROGS = variance

all : $(PROGS)

install : $(PROGS)
	ar rv $(libdir)/libds.a \
        variance_module.o
	cp $(PROGS) $(bindir)
	cp variance_module.h $(includedir)

variance : variance_module.o

variance_module.o : 

clean : 
	rm $(PROGS) \
	   variance_module.o
