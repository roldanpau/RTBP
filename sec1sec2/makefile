SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
#LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

PROGS = sec1sec2 sec1sec2_inv

all : $(PROGS)

install : $(PROGS) sec1sec2_module.o sec1sec2_module.h
	cp $(PROGS) $(bindir)
	ar rv $(libdir)/libds.a sec1sec2_module.o
	cp sec1sec2_module.h $(includedir)

sec1sec2: sec1sec2_module.o $(libdir)/libds.a

sec1sec2_inv: sec1sec2_module.o $(libdir)/libds.a

sec1sec2_module.o : $(includedir)/rtbp.h $(includedir)/hinv.h $(includedir)/prtbp2.h 

clean : 
	rm $(PROGS) \
	   sec1sec2_module.o \
