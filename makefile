# top-level Makefile
#
SHELL = /bin/sh

# Export these variables.
# Variables used by implicit rules like these are NOT normally passed down.
# To override variables used by the sub-make, one has to use the -e switch 
# to make, see below.
prefix = $(HOME)
exec_prefix = $(prefix)
export bindir = $(exec_prefix)/bin/rtbp
export includedir = $(prefix)/include/rtbp
export libdir = $(exec_prefix)/lib/rtbp
export LDFLAGS = -O3 -L$(HOME)/lib/rtbp
export CFLAGS = -O3 -DNDEBUG -I$(HOME)/include/rtbp

DIRS = rtbp taylor frtbp section hinv cardel prtbp_del_car prtbp utils \
       intersec_del_car prtbp_noloops errmfld invmfld invmfld_del_car \
       rtbp_del frtbp_red hinv_del frtbp_del prtbp_del \
       inner_circ outer_circ \
       initcond initcond_apo dprtbp portbp portbp_apo\
       sec1sec2 \
       hyper \
       approxint_del_car \
       inner_ell_stoch outer_ell_stoch \
	   approxint intersec splitting\
       trtbp

# the sets of directories to do various things in
BUILDDIRS = $(DIRS:%=build-%)
INSTALLDIRS = $(DIRS:%=install-%)
CLEANDIRS = $(DIRS:%=clean-%)
TESTDIRS = $(DIRS:%=test-%)

all: $(BUILDDIRS)
$(DIRS): $(BUILDDIRS)

# Configure makefile for debug/release builds
debug: export LDFLAGS = -g -L$(HOME)/lib/rtbp
debug: export CFLAGS = -g -I$(HOME)/include/rtbp
debug: $(BUILDDIRS)

# For each subdir, determine the subdir name (by stripping off the 
# "install-") and do a "make" in that dir.
#
# make "install" right after making "build" in each dir.
#
# To override variables CFLAGS/LDFLAGS used by the sub-make, one has to use the
# -e switch.
$(BUILDDIRS):
	$(MAKE) -e -C $(@:build-%=%) 
	$(MAKE) -e -C $(@:build-%=%) install

# build dependencies
build-taylor: install-rtbp
build-frtbp: install-taylor
build-prtbp: install-prtbp_noloops
build-cardel:install-hinv install-utils
build-prtbp_del_car: install-cardel install-section install-frtbp \
    install-prtbp_del install-utils
build-intersec_del_car: install-utils install-prtbp install-prtbp_del_car \
	install-errmfld
build-errmfld: install-prtbp_noloops
build-invmfld: install-errmfld
build-invmfld_del_car: install-errmfld install-invmfld \
	install-approxint_del_car
build-frtbp_red: install-rtbp_del
build-frtbp_del: install-rtbp_del
build-prtbp_del: install-frtbp_del install-hinv_del
build-inner_circ: install-frtbp_red
build-outer_circ: install-frtbp_del install-prtbp_del install-inner_circ
build-portbp: install-initcond install-dprtbp
build-portbp_apo: install-initcond_apo install-dprtbp
build-sec1sec2: install-prtbp

install: $(INSTALLDIRS)

$(INSTALLDIRS):
	$(MAKE) -e -C $(@:install-%=%) install

# install dependencies
install-rtbp : build-rtbp
install-taylor: build-taylor
install-frtbp: build-frtbp
install-section: build-section
install-hinv: build-hinv
install-cardel: build-cardel
install-prbtp_del_car: build-prtbp_del_car
install-prbtp: build-prtbp
install-utils: build-utils
install-intersec_del_car: build-intersec_del_car
install-prtbp_noloops: build-prtbp_noloops
install-errmfld: build-errmfld
install-invmfld: build-invmfld
install-invmfld_del_car: build-invmfld_del_car
install-rtbp_del : build-rtbp_del
install-frtbp_red : build-frtbp_red
install-hinv_del : build-hinv_del
install-frtbp_del : build-frtbp_del
install-prtbp_del : build-prtbp_del
install-inner_circ: build-inner_circ
install-outer_circ: build-outer_circ
install-initcond: build-initcond
install-initcond_apo: build-initcond_apo
install-dprtbp: build-dprtbp
install-portbp: build-portbp
install-portbp_apo: build-portbp_apo
install-sec1sec2: build-portbp
install-hyper: build-hyper
install-approxint_del_car: build-approxint_del_car
install-inner_ell_stoch: build-inner_ell_stoch
install-outer_ell_stoch: build-outer_ell_stoch
install-approxint: build-approxint
install-intersec: build-intersec
install-splitting: build-splitting
install-trtbp: build-trtbp

test: $(TESTDIRS) all
$(TESTDIRS): 
	$(MAKE) -C $(@:test-%=%) test

clean: $(CLEANDIRS) cleanlib

$(CLEANDIRS): 
	$(MAKE) -e -C $(@:clean-%=%) clean

cleanlib:
	rm $(libdir)/libds.a

.PHONY: $(DIRS)
.PHONY: $(BUILDDIRS)
.PHONY: $(INSTALLDIRS)
.PHONY: $(TESTDIRS)
.PHONY: $(CLEANDIRS)
.PHONY: all install clean cleanlib test

