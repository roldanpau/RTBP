SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
LDLIBS = -lm -lgsl -lgslcblas -lds

PROGS = hyperdel hyperdels

all : $(PROGS)

install : $(PROGS) hyperdel.o hyperdel.h
	ar rv $(libdir)/libds.a hyperdel.o 
	cp hyperdel.h $(includedir)
	cp $(PROGS) $(bindir)

hyperdels : hyperdels_main.o hyperdel.o $(libdir)/libds.a
	$(CC) -o hyperdels $(LDLIBS) $(CFLAGS) hyperdels_main.o hyperdel.o \
	$(libdir)/libds.a

hyperdels_main.o : 

hyperdel : hyperdel_main.o hyperdel.o $(libdir)/libds.a
	$(CC) -o hyperdel $(LDLIBS) $(CFLAGS) hyperdel_main.o hyperdel.o \
	$(libdir)/libds.a

hyperdel_main.o : 

hyperdel.o : $(includedir)/dprtbp_2d.h

clean : 
	rm hyperdels hyperdels_main.o \
	   hyperdel hyperdel_main.o hyperdel.o 
