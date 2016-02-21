SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
LDFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : dprtbp2 dprtbp2_test dprtbp2_2d dprtbp2_2d_test

install : dprtbp2.o dprtbp2.h dprtbp2 dprtbp2_2d.o dprtbp2_2d.h dprtbp2_2d
	ar rv $(libdir)/libds.a dprtbp2.o dprtbp2_2d.o
	cp dprtbp2.h dprtbp2_2d.h $(includedir)
	cp dprtbp2 dprtbp2_2d $(bindir)

dprtbp2 : dprtbp2_main.o dprtbp2.o $(libdir)/libds.a
	$(CC) -o dprtbp2 $(LDLIBS) $(CFLAGS) dprtbp2_main.o dprtbp2.o \
	$(libdir)/libds.a

dprtbp2_test : dprtbp2_test.o dprtbp2.o $(libdir)/libds.a
	$(CC) -o dprtbp2_test $(LDLIBS) $(CFLAGS) dprtbp2_test.o dprtbp2.o \
	$(libdir)/libds.a

dprtbp2_test.o : $(includedir)/rtbp.h $(includedir)/prtbp.h dprtbp2.h

dprtbp2_main.o : $(includedir)/rtbp.h

dprtbp2.o : $(includedir)/prtbp.h $(includedir)/frtbp.h $(includedir)/rtbp.h

dprtbp2_2d : dprtbp2_2d_main.o dprtbp2_2d.o dprtbp2.o $(libdir)/libds.a
	$(CC) -o dprtbp2_2d $(LDLIBS) $(CFLAGS) dprtbp2_2d_main.o dprtbp2_2d.o \
	dprtbp2.o $(libdir)/libds.a

dprtbp2_2d_main.o : dprtbp2_2d.h

dprtbp2_2d.o : $(includedir)/rtbp.h $(includedir)/prtbp.h \
   $(includedir)/hinv.h dprtbp2.h dprtbp2_2d.h

dprtbp2_2d_test : dprtbp2_2d_test.o dprtbp2_2d.o dprtbp2.o $(libdir)/libds.a
	$(CC) -o dprtbp2_2d_test $(LDLIBS) $(CFLAGS) dprtbp2_2d_test.o \
	dprtbp2_2d.o dprtbp2.o $(libdir)/libds.a

dprtbp2_2d_test.o : $(includedir)/prtbp_2d.h dprtbp2_2d.h

clean : 
	rm dprtbp2 dprtbp2_main.o dprtbp2.o \
	   dprtbp2_2d dprtbp2_2d_main.o dprtbp2_2d.o \
	   dprtbp2_test dprtbp2_test.o \
	   dprtbp2_2d_test dprtbp2_2d_test.o
