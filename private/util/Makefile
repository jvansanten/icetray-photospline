# definitions

CC = cc
AR = ar

LIBSRC = ../lib/bspline.c ../lib/fitstable.c
SRC = dumpsplinefits.c evalsplinefits.c
LIBS = -lm -L/usr/local/lib -lcfitsio

# -----

CFLAGS = -Wall -Werror -Wno-unused -march=native -O2 -I../lib -I/usr/local/include

RM = /bin/rm
.SUFFIXES: .o .c

LIBOBJ=$(LIBSRC:.c=.o)
OBJ=$(SRC:.c=.o)

dumpsplinefits: $(LIBOBJ) dumpsplinefits.o
	$(CC) -o dumpsplinefits $(LIBOBJ) dumpsplinefits.o $(LIBS)
evalsplinefits: $(LIBOBJ) evalsplinefits.o
	$(CC) -o evalsplinefits $(LIBOBJ) evalsplinefits.o $(LIBS)

$(LIBOBJ): %.o: %.c
	$(CC) -g -c $(CFLAGS) $< -o $@
$(OBJ): %.o: %.c
	$(CC) -g -c $(CFLAGS) $< -o $@

all: dumpsplinefits evalsplinefits

clean:
	$(RM) -rf *.o dumpsplinefits evalsplinefits
