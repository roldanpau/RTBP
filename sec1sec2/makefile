SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3 #-g -O
#LDFLAGS = -O3
LDLIBS = -lds -lm -lgsl -lgslcblas

PROGS = sec1sec2 sec1sec2_inv sec1sec2_2d sec1sec2_inv_2d

all : $(PROGS)

install : $(PROGS) \
    sec1sec2_module.o sec1sec2_module.h \
    sec1sec2_2d_module.o sec1sec2_2d_module.h 
	cp $(PROGS) $(bindir)
	ar rv $(libdir)/libds.a \
        sec1sec2_module.o \
        sec1sec2_2d_module.o
	cp sec1sec2_module.h sec1sec2_2d_module.h $(includedir)

sec1sec2: sec1sec2_module.o

sec1sec2_inv: sec1sec2_module.o

sec1sec2_module.o : $(includedir)/rtbp.h $(includedir)/hinv.h \
    $(includedir)/prtbp.h 

sec1sec2_2d: sec1sec2_2d_module.o

sec1sec2_inv_2d: sec1sec2_2d_module.o

sec1sec2_2d_module.o : $(includedir)/rtbp.h $(includedir)/hinv.h \
    $(includedir)/prtbp.h 

clean : 
	rm $(PROGS) \
	   sec1sec2_module.o \
	   sec1sec2_2d_module.o 
