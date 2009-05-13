#!/bin/sh

# This should probably be a makefile

CC_BOILERPLATE="splineutil.c ../lib/bspline.c -I/usr/local/include -L/usr/local/lib -lcholmod -lm -lgoto -lcolamd -lamd -I../lib -Wall"
CC="gcc"

echo 'Check the following against the python implementation:'
echo 'Testing kronecker product of rhotest x rhotest:'
$CC -o krontest krontest.c $CC_BOILERPLATE && time ./krontest < rhotest.m
echo 'Testing box of rhotest and rhotest:'
$CC -o boxtest boxtest.c $CC_BOILERPLATE && time ./boxtest < rhotest.m
echo 'Testing rho with the 3x3x3 array (0 1 2 ... 26) and rhotest along dimension 2:'
$CC -o rhotest rhotest.c $CC_BOILERPLATE && time ./rhotest < rhotest.m

