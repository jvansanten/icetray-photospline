# definitions

CC = cc
BLAS = goto
LAPACK = lapack
LIBDIRS = -L/usr/local/lib -L/usr/local/lib/gcc-4.3.4
INCDIRS = -I/usr/local/include -I/usr/local/include/python2.6

SRC = glam.c splineutil.c pyglam.c ../lib/bspline.c
LIBS = -lm -l$(BLAS) -lcholmod -lcolamd -lamd -lspqr -lstdc++ -llapack -lgfortran

# -----

CFLAGS = -g -fPIC -Wall -Werror -Wno-unused -march=native -I../lib $(INCDIRS)
LDFLAGS = -shared $(LIBDIRS)

RM = /bin/rm
.SUFFIXES: .o .c

OBJ=$(SRC:.c=.o)

spglam.so: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): %.o: %.c
	$(CC) -g -c $(CFLAGS) $< -o $@

all: spglam.so

clean:
	$(RM) -rf *.o *.so ../lib/*.o
