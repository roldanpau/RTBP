SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
#LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

PROGS = Lbound

all : $(PROGS)

install : $(PROGS)
	ar rv $(libdir)/libds.a \
        Lbound_module.o
	cp $(PROGS) $(bindir)
	cp Lbound_module.h $(includedir)

Lbound : Lbound_module.o

Lbound_module.o : $(includedir)/utils_module.h $(includedir)/frtbp.h \
$(includedir)/cardel.h

clean : 
	rm $(PROGS) \
	   Lbound_module.o
