SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm 

all : initcond

install : initcond
	cp initcond $(bindir)
	ar rv $(libdir)/libds.a initcond.o
	cp initcond.h $(includedir)

initcond : initcond_main.o initcond.o $(libdir)/libds.a
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

initcond_main.o : initcond.h

initcond.o :

clean : 
	rm initcond initcond_main.o initcond.o
