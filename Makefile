################################################################
# $Id: Makefile,v 1.2 2002/03/02 07:26:16 frolov Exp $
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

VERSION = 0.1



################################################################
# Configuration
################################################################

# Define if you want debugging code compiled in
#DEBUG = 1

# Define if you have FFTW library (for distortion analysis)
FFTW = 1

FFTINCS = -I$(HOME)/local/include
FFTLIBS = -L$(HOME)/local/lib -lrfftw -lfftw


# Compiler flags
ifdef DEBUG
CFLAGS	= -m486 -g -pg -Wall
DEFINES	= -DDEBUG -DVERSION='"$(VERSION) (debugging on)"' -DDATADIR='"$(DATADIR)"'
LDFLAGS	= -pg
else
CFLAGS	= -m486 -O6 -fomit-frame-pointer -funroll-loops
DEFINES	= -DVERSION='"$(VERSION)"' -DDATADIR='"$(DATADIR)"'
LDFLAGS	=
endif

INCLUDE	= 
LIBS	= -lm

ifdef FFTW
DEFINES:= -DFFTW $(DEFINES)
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


# Targets
.PHONY: all clean distclean install

all: $(bins)

clean:
	rm -f $(objs) gmon.out

distclean: clean
	rm -f $(bins)

install: $(bins)
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -s $(bins) $(BINDIR)
	$(INSTALL) curvecaptor $(BINDIR)


# Dependencies
tubefit: tubefit.o



################################################################
# Implicit rules
################################################################

# Production rules
%.o: %.c $(hdrs)
	$(CC) $(ALLFLAGS) -c $< -o $@

%: %.o
	$(CC) $(LDFLAGS) $^ $(LIBS) -o $@
