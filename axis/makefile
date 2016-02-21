SHELL = /bin/sh
prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib
CFLAGS = -O3
LDLIBS = -lm 

all : axis

install : axis
	#cp axis $(bindir)
	#ar rv $(libdir)/libds.a axis.o
	#cp axis.h $(includedir)

axis : axis_main.o 
	$(CC) -o axis $(LDLIBS) $(CFLAGS) axis_main.o

axis_main.o : 

clean : 
	rm axis axis_main.o
