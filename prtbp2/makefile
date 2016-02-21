SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

PROGS = prtbp2 prtbp2_inv prtbp2_2d prtbp2_2d_inv prtbp2s_2d

all : $(PROGS)

install : $(PROGS) \
   prtbp2.o prtbp2_2d.o \
   prtbp2.h prtbp2_2d.h
	cp $(PROGS) $(bindir)
	ar rv $(libdir)/libds.a prtbp2.o prtbp2_2d.o
	cp prtbp2.h prtbp2_2d.h $(includedir)

prtbp2 : prtbp2_main.o prtbp2.o $(libdir)/libds.a

prtbp2_main.o : $(includedir)/rtbp.h prtbp2.h

prtbp2.o : $(includedir)/frtbp.h $(includedir)/rtbp.h

prtbp2_inv : prtbp2_inv_main.o prtbp2.o $(libdir)/libds.a
	$(CC) -o prtbp2_inv $(LDLIBS) $(CFLAGS) prtbp2_inv_main.o prtbp2.o \
	   $(libdir)/libds.a

prtbp2_inv_main.o : $(includedir)/rtbp.h prtbp2.h

prtbp2_2d : prtbp2_2d_main.o prtbp2_2d.o prtbp2.o $(libdir)/libds.a
#	$(CC) -o prtbp2_2d $(LDLIBS) $(CFLAGS) prtbp2_2d_main.o prtbp2_2d.o \
#	prtbp2.o

prtbp2_2d_main.o : prtbp2_2d.h

prtbp2_2d.o : $(includedir)/rtbp.h $(includedir)/hinv.h prtbp2.h 

prtbp2_2d_inv : prtbp2_2d_inv_main.o prtbp2_2d.o prtbp2.o $(libdir)/libds.a
	$(CC) -o prtbp2_2d_inv $(LDLIBS) $(CFLAGS) prtbp2_2d_inv_main.o \
	   prtbp2_2d.o prtbp2.o $(libdir)/libds.a

prtbp2_2d_inv_main.o : prtbp2_2d.h

prtbp2s_2d : prtbp2s_2d_main.o prtbp2_2d.o prtbp2.o $(libdir)/libds.a
	$(CC) -o prtbp2s_2d $(LDLIBS) $(CFLAGS) prtbp2s_2d_main.o prtbp2_2d.o \
	prtbp2.o $(libdir)/libds.a

prtbp2s_2d_main.o : prtbp2_2d.h

clean : 
	rm $(PROGS) \
	   prtbp2_main.o prtbp2.o \
	   prtbp2_2d_main.o prtbp2_2d.o \
	   prtbp2s_2d_main.o \
	   prtbp2_inv_main.o \
	   prtbp2_2d_inv_main.o  \
