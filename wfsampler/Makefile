# definitions

CC = clang
LIBDIRS = -L/usr/local/lib
INCDIRS = -I/usr/local/include -I/usr/local/include/python2.6 -I/usr/local/lib/python2.6/site-packages/numpy/core/include

SRC = ../lib/splinepdf.c wfsample_numpy.c ../lib/bspline.c ../lib/fitstable.c
LIBS = -lm -lgsl -lcfitsio -lgslcblas

# -----

CFLAGS = -g -fPIC -Wall -Werror -Wno-unused -march=native -I../lib $(INCDIRS) -DUSE_CHOLMOD
LDFLAGS = -shared $(LIBDIRS)

RM = /bin/rm
.SUFFIXES: .o .c

OBJ=$(SRC:.c=.o)

wfsample.so: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): %.o: %.c
	$(CC) -g -c $(CFLAGS) $< -o $@

all: wfsample.so

clean:
	$(RM) -rf *.o *.so ../lib/*.o
