# Makefile.in for Easel's miniapps
# SVN $Id: Makefile.in 396 2009-10-06 18:09:50Z nawrockie $
SHELL = /bin/sh

## Installation targets
#
prefix      = /usr/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
bindir      = ${exec_prefix}/bin


## your compiler and compiler flags.
#
CC       = gcc -std=gnu99
CFLAGS   = -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing -msse2
DEFS     = -DHAVE_CONFIG_H
LDFLAGS  = 
LIBS     = -leasel   -lm

## where is Easel?
#
ESLLIBDIR = -L..
ESLINCDIR = -I..


## list of the miniapps to compile.
#
PROGS = esl-afetch\
	esl-alimanip\
	esl-alimap\
	esl-alistat\
	esl-compalign\
	esl-compstruct\
	esl-construct\
	esl-histplot\
	esl-reformat\
	esl-selectn\
	esl-seqstat\
	esl-sfetch\
	esl-shuffle\
	esl-ssdraw\
	esl-weight

## Targets defining how to make miniapp executables.
#
all: 	$(PROGS) 

${PROGS}: %: %.o  ../libeasel.a
	${CC} ${CFLAGS} ${DEFS} ${LDFLAGS} ${ESLLIBDIR} -o $@ $@.o ${OBJS} ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${DEFS} ${ESLINCDIR} -c $<

install:
	cp ${PROGS} ${DESTDIR}${bindir}/

clean:
	-rm ${PROGS}
	-rm *.o *~
	-rm TAGS

distclean:
	make clean
	-rm -f ${PROGS}
	-rm Makefile

TAGS:
	etags *.[ch] Makefile.in


################################################################
#    This copyrighted source code is freely distributed 
#    under the terms of the GNU General Public License. See
#    the files COPYRIGHT and LICENSE for details.
################################################################
