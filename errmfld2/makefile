SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : errmfld2.o

install : errmfld2.o
	#cp invmfld2 $(bindir)
	ar rv $(libdir)/libds.a errmfld2.o
	cp errmfld2.h $(includedir)

errmfld.o : $(includedir)/prtbp2_2d.h

clean : 
	rm errmfld2.o
