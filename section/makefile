SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : section.o

install : section.o section.h
	ar rv $(libdir)/libds.a section.o
	cp section.h $(includedir)

clean : 
	rm section.o
