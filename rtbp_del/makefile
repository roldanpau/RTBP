SHELL = /bin/sh

hostname := $(shell hostname)
ifeq ($(hostname),pau-PORTEGE-R700)
   prefix = /usr/local
   exec_prefix = $(prefix)
   bindir = $(exec_prefix)/bin
   includedir = $(prefix)/include
   libdir = $(exec_prefix)/lib
   CFLAGS = -O3 #-g
   LDLIBS = -lm -lgsl -lgslcblas #-lds
else
   prefix = $(HOME)
   exec_prefix = $(prefix)
   bindir = $(exec_prefix)/bin
   includedir = $(prefix)/include
   libdir = $(exec_prefix)/lib
   CFLAGS = -O3 -I$(GSL_INC) -I$(includedir) #-g
   LDFLAGS = -L$(GSL_LIB) -Wl,-rpath,$(GSL_LIB)
   LDLIBS = -lm -lgsl -lgslcblas #-lds
endif

all : rtbpdel

install : rtbpdel rtbpdel.o rtbpdel.h
	ar rv $(libdir)/libds.a rtbpdel.o
	cp rtbpdel.h $(includedir)
	cp rtbpdel $(bindir)

rtbpdel : rtbpdel_main.o rtbpdel.o
#	$(CC) -o frtbp $(LDLIBS) $(CFLAGS) frtbp_main.o frtbp.o rtbpdel.o

rtbpdel_main.o : rtbpdel.h

rtbpdel.o : rtbpdel.h

clean : 
	rm rtbpdel rtbpdel_main.o rtbpdel.o
