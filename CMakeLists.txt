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

SET_TARGET_PROPERTIES(photospline-evalsplinefits
        PROPERTIES
        COMPILE_FLAGS "-std=c99"
)

# We wants it fast, precious.
SET_TARGET_PROPERTIES(photospline
	PROPERTIES
	COMPILE_FLAGS "-O3 -fpeel-loops -fbranch-probabilities -ffast-math -std=c99"
	LINK_FLAGS "-fpeel-loops -fbranch-probabilities -ffast-math"
)
