#################################################################
# Infernal
# top level Makefile: installation, documentation
# CVS $Id: Makefile.in 2581 2008-09-18 13:37:59Z eddys $
#################################################################
#    This copyrighted source code is freely distributed 
#    under the terms of the GNU General Public License. See
#    the files COPYRIGHT and LICENSE for details.
#################################################################

# On most Unices, you can build the package without modifying anything 
#   by just typing "./configure; make".
#
# Post-configuration, you may still want to modify the following make variables:
#   BINDIR  - where the executables will be installed by a 'make install'
#   MANDIR  - where the man pages will be installed by a 'make install'
#   CC      - which compiler to use
#   CFLAGS  - compiler flags to use

# where you want things installed
# Sort of uses GNU coding standards. ${prefix} might be /usr/local.
# ${exec_prefix} gives you some flexibility for installing architecture
# dependent files (e.g. the programs): an example ${exec_prefix} might be
# /nfs/share/irix64/
#
prefix      = /usr/local
exec_prefix = ${prefix}
BINDIR      = ${exec_prefix}/bin
MANDIR      = ${datarootdir}/man
DATAROOTDIR = ${prefix}/share

# how to install the man pages; 
# cp is generally fine, unless you preformat your pages.
#
INSTMAN   = cp
MANSUFFIX = 1

# your compiler and compiler flags
#
CC        = gcc -std=gnu99
CFLAGS    = -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing
SIMDFLAGS = -msse2

#######
## You should not need to modify below this line.
## Some of it is concerned with maintenance of the development version
## and building the release (indeed, several commands will only work at
## Janelia)
#######
SHELL       = /bin/sh
PACKAGE     = Infernal
BASENAME    = infernal
RELEASE     = 1.0.2
RELEASEDATE = "October 2009"
COPYRIGHT   = "Copyright (C) 2009 HHMI Janelia Farm Research Campus"
LICENSELINE = "Freely distributed under the GNU General Public License (GPLv3)"
LICENSETAG  = gnu
COMPRESS    = gzip
REPOSITORY  = https://svn.janelia.org/eddylab/eddys/src/


# The program lists below for INFERNAL are not necessarily
# a complete manifest. They are the list of stable programs that the
# package will install. (The list in src/Makefile.in is the complete manifest
# of INFERNAL programs.)
#
PROGS = cmalign\
	cmbuild\
	cmcalibrate\
	cmemit\
	cmscore\
	cmsearch\
	cmstat

# Optional programs for rigorous filters
RFPROGS = cm2hmm cm2hmmsearch

# all: if we don't have precompiled binaries, compiles all packages, 
#      and moves the supported executables into binaries/ subdirectory.
#
#all: core $(RIGFILTERS)
all: core 

core:
	(cd easel;     make CC="$(CC)" CFLAGS="$(CFLAGS)"; make)
	(cd src;       make CC="$(CC)" CFLAGS="$(CFLAGS)"; make module)
	(cd testsuite; make CC="$(CC)" CFLAGS="$(CFLAGS)")

#.PHONY: $(RIGFILTERS)
#$(RIGFILTERS): core
#	if test -d $(RIGFILTERS); then (cd rigfilters; make); fi

# dist: build a new distribution directory in infernal-$RELEASE
#
dist:
	# Make sure we pointed the make dist to a subversion repository
	@if test x${REPOSITORY} = 'x'; then\
	  echo "REPOSITORY isn't set; use make REPOSITORY=... dist";\
	  exit 1;\
	fi
	# Delete old versions of the same release
	#
	@if test -d ../${BASENAME}-$(RELEASE);        then rm -rf ../${BASENAME}-$(RELEASE);        fi
	@if test -e ../${BASENAME}-$(RELEASE).tar;    then rm -f  ../${BASENAME}-$(RELEASE).tar;    fi
	@if test -e ../${BASENAME}-$(RELEASE).tar.Z;  then rm -f  ../${BASENAME}-$(RELEASE).tar.Z;  fi
	@if test -e ../${BASENAME}-$(RELEASE).tar.gz; then rm -f  ../${BASENAME}-$(RELEASE).tar.gz; fi
	# 
	# Export the release from Subversion control.
        # Because we use svn:externals, this automatically exports 
	# the ssdk, and Easel libraries as well.
	#
	svn export ${REPOSITORY} ../${BASENAME}-${RELEASE}
	# 
	# Make the configure scripts from configure.ac
	#
	(cd ../${BASENAME}-${RELEASE};        autoconf)
	(cd ../${BASENAME}-${RELEASE}/easel;  autoconf)
	#
	# We'll need a basic Makefiles for some steps; later
	# we'll delete them
	#
	(cd ../${BASENAME}-${RELEASE}/; ./configure)
	#
	# Make symlinks and set file permissions
	#
	(cd ../${BASENAME}-${RELEASE}/testsuite; ln -s ../easel/devkit/sqc .)
	find ../${BASENAME}-${RELEASE} -type f -exec chmod +r {} \;
	chmod +x ../${BASENAME}-${RELEASE}/easel/devkit/sqc
	chmod +x ../${BASENAME}-${RELEASE}/easel/devkit/autodoc
	#
	# Attach license stamps on all files that need 'em
	# (e.g. replace LICENSE keyword everywhere)
	#
	find ../${BASENAME}-${RELEASE} -type f -exec easel/devkit/sedition-pp LICENSE Licenses/$(LICENSETAG) {} \;
	#
	# Other keyword replacements
	#
	easel/devkit/sedition RELEASE ${RELEASE} RELEASEDATE ${RELEASEDATE}                        ../${BASENAME}-${RELEASE}/documentation/userguide/titlepage.tex
	easel/devkit/sedition RELEASE ${RELEASE} RELEASEDATE ${RELEASEDATE} COPYRIGHT ${COPYRIGHT} ../${BASENAME}-${RELEASE}/00README
	for manpage in ../${BASENAME}-${RELEASE}/documentation/manpages/*; do\
	   easel/devkit/sedition RELEASE ${RELEASE} RELEASEDATE ${RELEASEDATE} COPYRIGHT ${COPYRIGHT} PACKAGE ${PACKAGE} LICENSELINE ${LICENSELINE} $$manpage;\
	done
	#
	# Make the User's Guide, and copy to top level. (We'll
	# destroy the userguide source soon.)
	#
	(cd ../${BASENAME}-${RELEASE}/documentation/userguide; make pdf)
	cp ../${BASENAME}-${RELEASE}/documentation/userguide/main.pdf  ../${BASENAME}-${RELEASE}/Userguide.pdf
	# 
	# Finish the man pages off, adding necessary boilerplate.
	# Must be done *after* the user guide is created, since the user guide
	# steals and uses the pre-boilerplate man pages.
	# Include the optional rigfilters documentation too; whether the user
	# wants it or not is a configuration-time decision.
	#
	for prog in $(PROGS); do\
	   cat ../${BASENAME}-${RELEASE}/documentation/manpages/boilerplate.trailer >> ../${BASENAME}-${RELEASE}/documentation/manpages/$$prog.man;\
	done
	for prog in $(RFPROGS); do\
	   cat ../${BASENAME}-${RELEASE}/documentation/manpages/boilerplate.trailer >> ../${BASENAME}-${RELEASE}/documentation/manpages/$$prog.man;\
	done;
	-rm -f ../${BASENAME}-${RELEASE}/documentation/manpages/boilerplate.trailer
	#
	# Remove SVN-controlled files/directories that don't belong in the
	# distro
	#
	-rm -rf ../${BASENAME}-${RELEASE}/Licenses
	-rm -rf ../${BASENAME}-${RELEASE}/00CHECKLIST
	-rm -rf ../${BASENAME}-${RELEASE}/Bugs
	-rm -rf ../${BASENAME}-${RELEASE}/Manuscripts
	-rm -rf ../${BASENAME}-${RELEASE}/SUBPACKAGES
	-rm -rf ../${BASENAME}-${RELEASE}/00WISHLIST
	(cd ../${BASENAME}-${RELEASE}/easel; make implode)
	(cd ../${BASENAME}-${RELEASE}/;      make distclean)
	#
	# Packaging commands
	#
	(cd ..;	tar cvf ${BASENAME}-${RELEASE}.tar ${BASENAME}-${RELEASE})
	(cd ..; $(COMPRESS) ${BASENAME}-$(RELEASE).tar)

# "make check" runs a set of sqc exercises in the test suite, at level 4.
# To run faster, reduce the level (level 1 runs a few tests, level 4
# runs a lot). 
#
check:
	(cd testsuite; ./sqc -p ../easel/miniapps 4 exercises.sqc . ../src)

# "make devcheck" runs an additional set of development sqc exercises 
# in the test suite, at level 4. This will take a long time.
#
devcheck:
	(cd testsuite;\
	 ./sqc -p ../src/           4 esl-alimanip.sqc . ../easel/miniapps/;\
	 ./sqc -p ../easel/miniapps 4 cmalign.sqc      . ../src;\
	 ./sqc -p ../easel/miniapps 4 cmbuild.sqc      . ../src;\
	 ./sqc -p ../easel/miniapps 4 cmcalibrate.sqc  . ../src;\
	 ./sqc -p ../easel/miniapps 4 cmemit.sqc       . ../src;\
	 ./sqc -p ../easel/miniapps 4 cmscore.sqc      . ../src;\
	 ./sqc -p ../easel/miniapps 4 cmsearch.sqc     . ../src;\
	 ./sqc -p ../easel/miniapps 4 cmstat.sqc       . ../src;)

# "make install" installs the programs in BINDIR
#
install:
	mkdir -p ${BINDIR}
	for file in $(PROGS); do\
		cp src/$$file $(BINDIR)/;\
	done
#	if test -d $(RIGFILTERS); then\
#	   for file in $(RFPROGS); do\
#	      cp rigfilters/cm2hmm-1.0/$$file $(BINDIR)/;\
#	   done;\
#	fi

# "make clean" removes almost everything except configuration files
#              and binaries.
clean:
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out
	(cd testsuite;               make clean)
	(if cd documentation/userguide; then make clean; fi)
	(cd src;                     make clean)
	(cd easel;                   make clean)
#	(cd rigfilters;              make clean)

# "make distclean" leaves a pristine source distribution.
#
distclean:
	-rm -f  config.cache config.log config.status
	-rm -rf autom4te.cache
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out
	(cd testsuite;               make distclean)
	(if cd documentation/userguide; then make distclean; fi)
	(cd src;                     make distclean)
	(cd easel;                   make distclean)
#	(cd rigfilters;              make distclean)
	-rm Makefile

