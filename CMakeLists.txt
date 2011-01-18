#
# $Id: CMakeLists.txt 64301 2010-05-25 16:49:40Z jvansanten $
#

option(BUILD_PHOTOSPLINE_FITTER "Build libraries for reading and fitting Photonics tables." OFF)

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

i3_test_executable(test 
	private/test/*.cxx
	USE_TOOLS boost log4cplus
	USE_PROJECTS photospline
)

if(BUILD_PHOTOSPLINE_FITTER)
	i3_add_library(photo2numpy
		private/photo2numpy/photo2numpy.c
		USE_TOOLS photonics python numpy
		INSTALL_DESTINATION lib/icecube/photospline
		NOT_INSPECTABLE
		MODULE
	)
	add_custom_command(TARGET photo2numpy
		PRE_LINK
		COMMAND mkdir -p ${CMAKE_BINARY_DIR}/lib/icecube/photospline
	)
	set_target_properties(photo2numpy PROPERTIES
		PREFIX ""
		OUTPUT_NAME photo2numpy
		LIBRARY_OUTPUT_DIRECTORY
		${LIBRARY_OUTPUT_PATH}/icecube/photospline
	)

	i3_add_library(spglam
		private/cfitter/cholesky_solve.c private/cfitter/glam.c
		private/cfitter/nnls.c private/cfitter/pyglam.c
		private/cfitter/splineutil.c
		USE_TOOLS blas python numpy suitesparse
		USE_PROJECTS photospline
		INSTALL_DESTINATION lib/icecube/photospline
		NOT_INSPECTABLE
		MODULE
	)
	add_custom_command(TARGET spglam
		PRE_LINK
		COMMAND mkdir -p ${CMAKE_BINARY_DIR}/lib/icecube/photospline
	)
	set_target_properties(spglam PROPERTIES
		PREFIX ""
		OUTPUT_NAME spglam 
		LIBRARY_OUTPUT_DIRECTORY
		${LIBRARY_OUTPUT_PATH}/icecube/photospline
	)
	add_custom_command(TARGET spglam
		PRE_LINK
		COMMAND mkdir -p ${CMAKE_BINARY_DIR}/lib/icecube/photospline
	)
endif(BUILD_PHOTOSPLINE_FITTER)
