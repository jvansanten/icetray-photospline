# definitions

CC = cc
#INCDIRS = -I/usr/local/include -I/usr/local/include/python2.6
INCDIRS = -I/afs/ifh.de/group/amanda/software/RHEL_5.0_amd64/include/python2.5 -I/afs/ifh.de/group/amanda/software/RHEL_5.0_amd64/lib/python2.5/site-packages/numpy/core/include -I/afs/ifh.de/user/n/nwhiteho/photonics/lib -I/afs/ifh.de/user/n/nwhiteho/photonics/level2
PHOTONICS_LIB = /afs/ifh.de/user/n/nwhiteho/photonics/lib/.libs/libphotonics.a \
	/afs/ifh.de/user/n/nwhiteho/photonics/level2/.libs/liblevel2amasim.a

SRC = photo2numpy.c
LIBS = -lm $(PHOTONICS_LIB)

# -----

CFLAGS = -O3 -fPIC -Wall -Werror -I../lib $(INCDIRS)
LDFLAGS = -shared 

RM = /bin/rm
.SUFFIXES: .o .c

OBJ=$(SRC:.c=.o)

photo2numpy.so: $(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

$(OBJ): %.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

all: photo2numpy.so

clean:
	$(RM) -rf *.o *.so ../lib/*.o
