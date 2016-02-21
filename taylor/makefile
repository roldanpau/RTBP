CC=gcc
CFLAGS=-O3 -Wall -pedantic
LFLAGS=-lm -s

TAYLOR=taylor

prefix = /usr/local
exec_prefix = $(prefix)
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib



all: prova prova3d

prova: prova.o taylor2dv.o
	$(CC) $(CFLAGS) prova.o taylor2dv.o -o prova $(LFLAGS)

prova.o: prova.c taylor2dv.h
	$(CC) $(CFLAGS) -c prova.c

prova3d.o: prova3d.c taylor3dv.h
	$(CC) $(CFLAGS) -c prova3d.c

taylor2dv.o: taylor2dv.c taylor2dv.h
	$(CC) $(CFLAGS) -c taylor2dv.c

taylor2dv.c: rtbp2dv.in
	$(TAYLOR) -name rtbp2dv -headername taylor2dv.h -o taylor2dv.c -jet \
	   -step -sqrt rtbp2dv.in

prova3d: prova3d.o taylor3dv.o
	$(CC) $(CFLAGS) prova3d.o taylor3dv.o -o prova3d $(LFLAGS)

taylor3dv.o: taylor3dv.c taylor3dv.h
	$(CC) $(CFLAGS) -c taylor3dv.c

taylor3dv.c: rtbp3dv.in
	$(TAYLOR) -name rtbp3dv -headername taylor3dv.h -o taylor3dv.c -jet \
	   -step -sqrt rtbp3dv.in

taylor2dv.h:
	$(TAYLOR) -name rtbp2dv -header -o taylor2dv.h rtbp2dv.in

taylor3dv.h:
	$(TAYLOR) -name rtbp3dv -header -o taylor3dv.h rtbp3dv.in

install: taylor2dv.o taylor2dv.h
	ar rv $(libdir)/libds.a taylor2dv.o
	cp taylor2dv.h $(includedir)


clean:
	rm -f prova prova3d
	rm -f *.o
	rm -f taylor2dv.h taylor3dv.h taylor2dv.c taylor3dv.c
