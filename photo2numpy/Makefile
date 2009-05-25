# definitions

CC = cc
#INCDIRS = -I/usr/local/include -I/usr/local/include/python2.6
INCDIRS = -I/usr/include/python2.3 -I/home/nwhitehorn/i3-software/x86_64/python/numpy/core/include -I/home/nwhitehorn/photonics/lib
PHOTONICS_LIB = /home/nwhitehorn/photonics/lib/.libs/libphotonics.a

SRC = photo2numpy.c
LIBS = -lm $(PHOTONICS_LIB)

# -----

CFLAGS = -O3 -fPIC -Wall -Werror -march=opteron -I../lib $(INCDIRS)
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
