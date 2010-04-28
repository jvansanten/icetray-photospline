# definitions

CC = clang
BLAS = goto2
LAPACK = lapack
LIBDIRS = -L/usr/local/lib -L/usr/local/lib/gcc44
#INCDIRS = -I/usr/local/include -I/usr/local/include/python2.6
INCDIRS = -I/usr/local/include -I/usr/local/include/python2.6 -I/usr/local/lib/python2.6/site-packages/numpy/core/include

SRC = glam.c nnls.c splineutil.c pyglam.c ../lib/bspline.c
LIBS = -lm -l$(BLAS) -lcamd -lcolamd -lamd -lstdc++ -llapack -lccolamd -lmetis -lgfortran -lcholmod -lspqr

# -----

CFLAGS = -g -fPIC -Wall -Werror -Wno-unused -march=native -I../lib $(INCDIRS) -DUSE_CHOLMOD
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
