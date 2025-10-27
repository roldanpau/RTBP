SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
#LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

PROGS = ebound

all : $(PROGS)

install : $(PROGS)
	ar rv $(libdir)/libds.a \
        ebound_module.o
	cp $(PROGS) $(bindir)
	cp ebound_module.h $(includedir)

ebound : ebound_module.o

ebound_module.o : $(includedir)/utils_module.h $(includedir)/frtbp.h \
$(includedir)/cardel.h

clean : 
	rm $(PROGS) \
	   ebound_module.o
