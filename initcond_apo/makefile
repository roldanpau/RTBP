SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDFLAGS = -O3
LDLIBS = -lm 

all : initcond_apo

install : initcond_apo
	cp initcond_apo $(bindir)
	ar rv $(libdir)/libds.a initcond_apo.o
	cp initcond_apo.h $(includedir)

initcond_apo : initcond_apo_main.o initcond_apo.o
#	$(CC) -o prtbp $(LDLIBS) $(CFLAGS) prtbp_main.o prtbp.o

initcond_apo_main.o : initcond_apo.h

initcond_apo.o :

clean : 
	rm initcond_apo initcond_apo_main.o initcond_apo.o
