SHELL = /bin/sh

hostname := $(shell hostname)
ifeq ($(hostname),deepthought2)
   prefix = $(HOME)
   exec_prefix = $(prefix)
   bindir = $(exec_prefix)/bin
   includedir = $(prefix)/include
   libdir = $(exec_prefix)/lib
   CFLAGS = -O3 -I$(GSL_INC)
   LDFLAGS = -L$(libdir)
   LDLIBS = -lm
else
   prefix = /usr/local
   exec_prefix = $(prefix)
   bindir = $(exec_prefix)/bin
   includedir = $(prefix)/include
   libdir = $(exec_prefix)/lib
   CFLAGS = -O3
   LDLIBS = -lm
endif

all : rtbp

install : rtbp rtbp.o rtbp.h
	ar rv $(libdir)/libds.a rtbp.o
	cp rtbp.h $(includedir)
	cp rtbp $(bindir)

rtbp : rtbp_main.o rtbp.o
#	$(CC) -o frtbp $(LDLIBS) $(CFLAGS) frtbp_main.o frtbp.o rtbp.o

rtbp_main.o : rtbp.h

rtbp.o : rtbp.h

clean : 
	rm rtbp rtbp_main.o rtbp.o
