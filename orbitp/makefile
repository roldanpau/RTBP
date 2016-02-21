SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm -lgsl -lgslcblas -lds

all : orbitp2 orbitp2_bwd

install : orbitp2 orbitp2_bwd orbitp2.o orbitp2.h
	cp orbitp2 orbitp2_bwd $(bindir)
	ar rv $(libdir)/libds.a orbitp2.o
	cp orbitp2.h $(includedir)

orbitp2 : orbitp2_main.o orbitp2.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

orbitp2_bwd : orbitp2_bwd_main.o orbitp2.o $(libdir)/libds.a
	$(CC) -o orbitp2_bwd $(LDLIBS) $(CFLAGS) orbitp2_bwd_main.o orbitp2.o \
	$(libdir)/libds.a

orbitp2_main.o : $(includedir)/prtbp2_2d.h orbitp2.h

orbitp2_bwd_main.o : $(includedir)/prtbp2_2d.h orbitp2.h

orbitp2.o : 

clean : 
	rm orbitp2 orbitp2_bwd orbitp2_main.o orbitp2_bwd_main.o orbitp2.o
