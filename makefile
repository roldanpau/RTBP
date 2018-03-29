# top-level Makefile
#
SHELL = /bin/sh

DIRS = rtbp taylor frtbp section hinv cardel prtbp_del_car prtbp utils \
       intersec_del_car prtbp_noloops errmfld invmfld invmfld_del_car \
       rtbp_del frtbp_red hinv_del frtbp_del prtbp_del \
       inner_circ outer_circ

# the sets of directories to do various things in
BUILDDIRS = $(DIRS:%=build-%)
INSTALLDIRS = $(DIRS:%=install-%)
CLEANDIRS = $(DIRS:%=clean-%)
TESTDIRS = $(DIRS:%=test-%)

all: $(BUILDDIRS)
$(DIRS): $(BUILDDIRS)

# For each subdir, determine the subdir name (by stripping off the 
# "install-") and do a "make" in that dir.
#
# make "install" right after making "build" in each dir.
$(BUILDDIRS):
	$(MAKE) -C $(@:build-%=%)
	sudo $(MAKE) -C $(@:build-%=%) install

# build dependencies
build-taylor: install-rtbp
build-frtbp: install-taylor
build-cardel:install-hinv
build-prtbp_del_car: install-cardel install-section install-frtbp
build-intersec_del_car: install-utils install-prtbp install-prtbp_del_car
build-errmfld: install-prtbp_noloops
build-invmfld: install-errmfld
build-invmfld_del_car: install-errmfld install-invmfld
build-frtbp_red: install-rtbp_del
build-prtbp_del: install-frtbp_del install-hinv_del
build-inner_circ: install-frtbp_red
build-outer_circ: install-frtbp_del install-prtbp_del install-inner_circ

install: $(INSTALLDIRS)

$(INSTALLDIRS):
	sudo $(MAKE) -C $(@:install-%=%) install

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
install-prbtp_noloops: build-prtbp_noloops
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

test: $(TESTDIRS) all
$(TESTDIRS): 
	$(MAKE) -C $(@:test-%=%) test

clean: $(CLEANDIRS) cleanlib

$(CLEANDIRS): 
	$(MAKE) -C $(@:clean-%=%) clean

cleanlib:
	sudo rm /usr/local/lib/libds.a

.PHONY: $(DIRS)
.PHONY: $(BUILDDIRS)
.PHONY: $(INSTALLDIRS)
.PHONY: $(TESTDIRS)
.PHONY: $(CLEANDIRS)
.PHONY: all install clean cleanlib test
