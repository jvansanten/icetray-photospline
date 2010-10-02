# definitions

CC = cc
AR = ar

LIBSRC = ../lib/bspline.c ../lib/fitstable.c ../lib/bspline_multi.c
SRC = dumpsplinefits.c evalsplinefits.c
LIBS = -lm -L/usr/local/lib -L${I3_PORTS}/lib -lcfitsio

# -----

CFLAGS = -ffast-math -fpeel-loops -fbranch-probabilities -Wall -Werror -Wno-unused -march=native -O3 -I../lib -I/usr/local/include -I../../public -I${SROOT}/include -I${I3_PORTS}/include -std=c99
#CFLAGS = -fprofile-generate -ffast-math -fpeel-loops -fbranch-probabilities -Wall -Werror -Wno-unused -march=native -O3 -I../lib -I/usr/local/include -I../../public -I${SROOT}/include -I${I3_PORTS}/include -std=c99

RM = /bin/rm
.SUFFIXES: .o .c

LIBOBJ=$(LIBSRC:.c=.o)
OBJ=$(SRC:.c=.o)

dumpsplinefits: $(LIBOBJ) dumpsplinefits.o
	$(CC) $(CFLAGS) -o dumpsplinefits $(LIBOBJ) dumpsplinefits.o $(LIBS)
evalsplinefits: $(LIBOBJ) evalsplinefits.o
	$(CC) $(CFLAGS) -o evalsplinefits $(LIBOBJ) evalsplinefits.o $(LIBS)

$(LIBOBJ): %.o: %.c
	$(CC) -g -c $(CFLAGS) $< -o $@
$(OBJ): %.o: %.c
	$(CC) -g -c $(CFLAGS) $< -o $@

all: dumpsplinefits evalsplinefits

clean:
	$(RM) -rf *.o dumpsplinefits evalsplinefits
