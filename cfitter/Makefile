# definitions

CC = cc
BLAS = goto
LAPACK = lapack
LIBDIRS = -L/usr/lib64 -L/afs/ifh.de/group/amanda/software/RHEL_5.0_amd64/lib -L/afs/ifh.de/user/n/nwhiteho/SuiteSparse/lib -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -rpath /afs/ifh.de/user/n/nwhiteho/SuiteSparse/lib
#INCDIRS = -I/usr/local/include -I/usr/local/include/python2.6
#INCDIRS = -I/usr/local/include -I/usr/local/include/python2.6 -I/usr/local/lib/python2.6/site-packages/numpy/core/include
INCDIRS = -I/afs/ifh.de/group/amanda/software/RHEL_5.0_amd64/include/python2.5 -I/afs/ifh.de/group/amanda/software/RHEL_5.0_amd64/lib/python2.5/site-packages/numpy/core/include -I/afs/ifh.de/user/n/nwhiteho/SuiteSparse

SRC = glam.c splineutil.c pyglam.c ../lib/bspline.c
LIBS = -lm -l$(BLAS) -lcholmod -lccolamd -lcolamd -lamd -lspqr -llapack -lgfortran -lstdc++ -lgfortranbegin /usr/lib64/libg2c.so.0

# -----

CFLAGS = -g -fPIC -Wall -Werror -Wno-unused -I../lib $(INCDIRS)
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
