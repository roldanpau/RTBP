SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
#LDFLAGS = -O3
LDLIBS = -lds -lm

PROGS = utils

all : $(PROGS) utils_module.o

utils : utils_module.o

install : utils_module.o utils_module.h
	ar rv $(libdir)/libds.a utils_module.o
	cp utils_module.h $(includedir)

clean : 
	rm $(PROGS) utils_module.o
