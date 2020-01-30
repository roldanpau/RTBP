CC=gcc
CFLAGS=-O3 -Wall -pedantic
LFLAGS=-lm -s

TAYLOR=taylor

prefix = $(HOME)
exec_prefix = $(prefix)
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib



all: prova2d prova2dv prova3dv

prova2d: prova2d.o taylor2d.o
	$(CC) $(CFLAGS) prova2d.o taylor2d.o -o prova2d $(LFLAGS)

prova2d.o: prova2d.c taylor2d.h
	$(CC) $(CFLAGS) -c prova2d.c

prova2dv: prova2dv.o taylor2dv.o
	$(CC) $(CFLAGS) prova2dv.o taylor2dv.o -o prova2dv $(LFLAGS)

prova2dv.o: prova2dv.c taylor2dv.h
	$(CC) $(CFLAGS) -c prova2dv.c

prova3dv: prova3dv.o taylor3dv.o
	$(CC) $(CFLAGS) prova3dv.o taylor3dv.o -o prova3dv $(LFLAGS)

prova3dv.o: prova3dv.c taylor3dv.h
	$(CC) $(CFLAGS) -c prova3dv.c

taylor2d.o: taylor2d.c taylor2d.h
	$(CC) $(CFLAGS) -c taylor2d.c

taylor2d.c: rtbp2d.in
	$(TAYLOR) -name rtbp2d -headername taylor2d.h -o taylor2d.c -jet \
	   -step -sqrt rtbp2d.in

taylor2dv.o: taylor2dv.c taylor2dv.h
	$(CC) $(CFLAGS) -c taylor2dv.c

taylor2dv.c: rtbp2dv.in
	$(TAYLOR) -name rtbp2dv -headername taylor2dv.h -o taylor2dv.c -jet \
	   -step -sqrt rtbp2dv.in

taylor3dv.o: taylor3dv.c taylor3dv.h
	$(CC) $(CFLAGS) -c taylor3dv.c

taylor3dv.c: rtbp3dv.in
	$(TAYLOR) -name rtbp3dv -headername taylor3dv.h -o taylor3dv.c -jet \
	   -step -sqrt rtbp3dv.in

taylor2d.h:
	$(TAYLOR) -name rtbp2d -header -o taylor2d.h rtbp2d.in

taylor2dv.h:
	$(TAYLOR) -name rtbp2dv -header -o taylor2dv.h rtbp2dv.in

taylor3dv.h:
	$(TAYLOR) -name rtbp3dv -header -o taylor3dv.h rtbp3dv.in

install: taylor2d.o taylor2d.h taylor2dv.o taylor2dv.h
	ar rv $(libdir)/libds.a taylor2d.o taylor2dv.o
	cp taylor2d.h taylor2dv.h $(includedir)


clean:
	rm -f prova2d prova2dv prova3dv
	rm -f *.o
	rm -f taylor2d.h taylor2dv.h taylor3dv.h taylor2d.c taylor2dv.c taylor3dv.c
