#
# $Id: CMakeLists.txt 64301 2010-05-25 16:49:40Z jvansanten $
#

i3_project(photospline
	PYTHON_DIR python
)

i3_add_library(photospline
	private/lib/*.c
	USE_TOOLS cfitsio gsl python
)

i3_executable(evalsplinefits
	private/util/evalsplinefits.c
	USE_TOOLS cfitsio
	USE_PROJECTS photospline)

i3_executable(dumpsplinefits
	private/util/dumpsplinefits.c
	USE_TOOLS cfitsio
	USE_PROJECTS photospline)

SET_TARGET_PROPERTIES(photospline-evalsplinefits
        PROPERTIES
        COMPILE_FLAGS "-std=c99"
)

include(vector_units.cmake)

# add flags to enable the vector units present on the host processor
# we'd love to use -march=native, but that's only available on GCC >= 4.2
if(NOT DEFINED PHOTOSPLINE_CFLAGS)
	SET(PHOTOSPLINE_CFLAGS "-O3 -fpeel-loops -fbranch-probabilities -ffast-math -std=c99")
	add_vector_flags(${CMAKE_C_COMPILER} PHOTOSPLINE_CFLAGS _junk)
	SET(PHOTOSPLINE_CFLAGS ${PHOTOSPLINE_CFLAGS} CACHE STRING "Compiler flags to enable vector units present on the host processor")
endif(NOT DEFINED PHOTOSPLINE_CFLAGS)

# We wants it fast, precious.
SET_TARGET_PROPERTIES(photospline
	PROPERTIES
	COMPILE_FLAGS ${PHOTOSPLINE_CFLAGS}
	LINK_FLAGS "-fpeel-loops -fbranch-probabilities -ffast-math"
)
