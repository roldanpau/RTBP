SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lgsl -lgslcblas -lds -lm

all : errmfld.o

install : errmfld.o
	ar rv $(libdir)/libds.a errmfld.o
	cp errmfld.h $(includedir)

errmfld.o : $(includedir)/prtbp_2d.h

clean : 
	rm errmfld.o
