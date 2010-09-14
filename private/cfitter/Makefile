# definitions

CC = gcc
BLAS = goto2
LAPACK = lapack
LIBDIRS = -L${SROOT}/lib
INCDIRS = -I${SROOT}/include -I../../public -I${SROOT}/lib/python2.6/site-packages/numpy/core/include -I${SROOT}/include/python2.6

SRC = glam.c cholesky_solve.c nnls.c splineutil.c pyglam.c ../lib/bspline.c
LIBS = -lm -l$(BLAS) -lcamd -lcolamd -lamd -lccolamd -lmetis -lgfortran -lcholmod

# -----

CFLAGS = -g -fPIC -Wall -Werror -Wno-unused -I../lib $(INCDIRS) -DUSE_CHOLMOD
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
