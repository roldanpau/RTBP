SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : errmflddel.o

install : errmflddel.o
	#cp invmfld2 $(bindir)
	ar rv $(libdir)/libds.a errmflddel.o
	cp errmflddel.h $(includedir)

errmfld.o : $(includedir)/prtbpdel_2d.h

clean : 
	rm errmflddel.o
