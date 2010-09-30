#
# $Id: CMakeLists.txt 64301 2010-05-25 16:49:40Z jvansanten $
#

i3_project(photospline)

i3_add_library(photospline
	private/lib/*.c
	USE_TOOLS cfitsio gsl python
)

i3_executable(evalsplinefits
	private/util/evalsplinefits.c
	USE_TOOLS cfitsio
	USE_PROJECTS photospline)

# We wants it fast, precious.
SET_TARGET_PROPERTIES(photospline
	PROPERTIES
	COMPILE_FLAGS "-O3 -ffast-math -std=c99"
)
