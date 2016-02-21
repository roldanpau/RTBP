SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : dprtbpdel_2d

install : dprtbpdel_2d.o dprtbpdel_2d.h dprtbpdel_2d
	ar rv $(libdir)/libds.a dprtbpdel_2d.o
	cp dprtbpdel_2d.h $(includedir)
	cp dprtbpdel_2d $(bindir)

dprtbpdel_2d : dprtbpdel_2d_main.o dprtbpdel_2d.o $(libdir)/libds.a
	$(CC) -o dprtbpdel_2d $(LDLIBS) $(CFLAGS) dprtbpdel_2d_main.o dprtbpdel_2d.o \
	$(libdir)/libds.a

dprtbpdel_2d_main.o : dprtbpdel_2d.h

dprtbpdel_2d.o : dprtbpdel_2d.h

clean : 
	rm dprtbpdel_2d dprtbpdel_2d_main.o dprtbpdel_2d.o
