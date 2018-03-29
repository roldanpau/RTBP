SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g
LDFLAGS = -O3 #-g

# PRG (2018/03/29): Moving -lm in front of -lds gives me trouble?!
LDLIBS = -lds -lgsl -lgslcblas -lm

all : inner_circ

install : inner_circ
	ar rv $(libdir)/libds.a inner_circ.o
	cp inner_circ.h $(includedir)
	cp inner_circ $(bindir)

#inner_circ : inner_circ_main.o inner_circ.o $(libdir)/libds.a
inner_circ : inner_circ_main.o inner_circ.o
#	$(CC) -o inner_circ $(LDLIBS) $(CFLAGS) inner_circ_main.o inner_circ.o rtbp.o

inner_circ_main.o : inner_circ.h

inner_circ.o : $(includedir)/rtbpdel.h $(includedir)/frtbpred.h

clean : 
	rm inner_circ \
	   inner_circ_main.o inner_circ.o
