SHELL = /bin/sh
prefix = $(HOME)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
LDLIBS = -lds -lm -lgsl -lgslcblas

all : trtbp

install : trtbp
	cp trtbp $(bindir)

trtbp : trtbp.o
#	cc -o trtbp -g -lm -lgsl -lgslcblas trtbp.o rtbp.o

#frtbp : frtbp_main.o frtbp.o libds.a(rtbp.o)
#	$(CC) -o frtbp $(LDLIBS) $(CFLAGS) frtbp_main.o frtbp.o rtbp.o

trtbp.o : $(includedir)/frtbp.h

clean : 
	rm trtbp trtbp.o
