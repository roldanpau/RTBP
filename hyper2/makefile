SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
LDLIBS = -lm -lgsl -lgslcblas -lds

RESULTS = hyper2 hyper2s

all : $(RESULTS)

install : $(RESULTS) hyper2.o hyper2.h
	cp $(RESULTS) $(bindir)
	ar rv $(libdir)/libds.a hyper2.o 
	cp hyper2.h $(includedir)

hyper2s : hyper2s.o hyper2.o $(libdir)/libds.a

hyper2s.o : 

hyper2 : hyper2_main.o hyper2.o $(libdir)/libds.a

hyper2_main.o : 

hyper2.o : $(includedir)/dprtbp_2d.h

clean : 
	rm $(RESULTS) \
	   hyper2s.o \
	   hyper2_main.o hyper2.o 
