################################################################
# $Id: Makefile,v 1.6 2005/06/02 22:14:48 afrolov Exp $
# Makefile for Curve Captor source code
################################################################

SHELL = /bin/sh

PREFIX  = /usr/local

BINDIR  = $(PREFIX)/bin
LIBDIR  = $(PREFIX)/lib
DATADIR = $(PREFIX)/share/valves


CC	= gcc
LD	= ld
AR	= ar -r
RANLIB	= ranlib
INSTALL = install

VERSION = 0.9.1



################################################################
# Configuration
################################################################

# Define if you want debugging code compiled in
# DEBUG = 1

# Define if you have getline() available
HAVE_GETLINE = 1

# Define if you have FFTW library (for distortion analysis)
# HAVE_FFTW = 1

FFTINCS = -I$(HOME)/local/include
FFTLIBS = -L$(HOME)/local/lib -lrfftw -lfftw


# Compiler flags
ifdef DEBUG
CFLAGS	= -g -pg -Wall
DEFINES	= -DDEBUG -DVERSION='"$(VERSION) (debugging on)"' -DDATADIR='"$(DATADIR)"'
LDFLAGS	= -pg
else
CFLAGS	= -O3 -march=i386 -mcpu=pentium4
DEFINES	= -DVERSION='"$(VERSION)"' -DDATADIR='"$(DATADIR)"'
LDFLAGS	=
endif

INCLUDE	= 
LIBS	= -lm

ifdef HAVE_GETLINE
DEFINES:= -DHAVE_GETLINE $(DEFINES)
endif

ifdef HAVE_FFTW
DEFINES:= -DHAVE_FFTW $(DEFINES)
INCLUDE:= $(FFTINCS) $(INCLUDE)
LIBS   := $(FFTLIBS) $(LIBS)
endif

ALLFLAGS = $(CFLAGS) $(DEFINES) $(INCLUDE)



################################################################
# Dependencies
################################################################

# Files
hdrs = $(wildcard *.h)
objs = $(wildcard *.o *.a)
bins = tubefit

dist = curvecaptor-$(VERSION).tar.gz


# Targets
.PHONY: all clean distclean install tarball

all: $(bins)

clean:
	rm -f $(objs) gmon.out

distclean: clean
	rm -f $(bins)

install: $(bins)
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -s $(bins) $(BINDIR)
	$(INSTALL) curvecaptor $(BINDIR)

tarball: $(dist)


# Dependencies
tubefit: tubefit.o

$(dist): README Makefile curvecaptor tubefit.c models.m4
	 tar cvf $(basename $@) $^; gzip -f -9 $(basename $@)



################################################################
# Implicit rules
################################################################

# Production rules
%.o: %.c $(hdrs)
	$(CC) $(ALLFLAGS) -c $< -o $@

%: %.o
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@
