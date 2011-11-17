#
# $Id: CMakeLists.txt 64301 2010-05-25 16:49:40Z jvansanten $
#

option(BUILD_PHOTOSPLINE_FITTER "Build libraries for reading and fitting Photonics tables." OFF)

i3_project(photospline
	PYTHON_DIR python
)

i3_add_library(photospline
	private/lib/bspline.c
	private/lib/bspline_multi.c
	private/lib/convolve.c
	private/lib/fitstable.c
	private/lib/splinepdf.c
	private/lib/I3SplineTable.cxx
	USE_TOOLS cfitsio gsl python
)

add_subdirectory(private/pybindings)

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

if(NOT DEFINED PHOTOSPLINE_CFLAGS)
	SET(PHOTOSPLINE_CFLAGS "-O3 -fpeel-loops -fbranch-probabilities -ffast-math -std=c99")
	# add flags to enable the vector units present on the host processor
	# we'd love to use -march=native, but that's only available on
	# GCC >= 4.2 and only on some architectures.
	if(CMAKE_SYSTEM_PROCESSOR MATCHES "^i.86$")
		SET(PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -msse2")
	elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(powerpc|ppc)")
		SET(PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -maltivec")
	elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^sparc")
		SET(PHOTOSPLINE_CFLAGS "${PHOTOSPLINE_CFLAGS} -mvis")
	endif(CMAKE_SYSTEM_PROCESSOR MATCHES "^i.86$")
	# Note: SSE2 is a requirement of amd64, so we don't worry about that
	SET(PHOTOSPLINE_CFLAGS ${PHOTOSPLINE_CFLAGS} CACHE STRING "Compiler flags specific to photospline")
endif(NOT DEFINED PHOTOSPLINE_CFLAGS)

# We wants it fast, precious.
SET_SOURCE_FILES_PROPERTIES(
	private/lib/bspline.c
	private/lib/bspline_multi.c
	private/lib/convolve.c
	private/lib/fitstable.c
	private/lib/splinepdf.c
	PROPERTIES
	COMPILE_FLAGS ${PHOTOSPLINE_CFLAGS}
)
SET_TARGET_PROPERTIES(photospline
	PROPERTIES
	LINK_FLAGS "-fpeel-loops -fbranch-probabilities -ffast-math"
)

i3_test_executable(test 
	private/test/*.cxx
	USE_TOOLS boost log4cplus
	USE_PROJECTS photospline
)

include(CheckIncludeFile)
CHECK_INCLUDE_FILE(level2_reader.h HAVE_LEVEL2_READER_H)

if (NOT HAVE_LEVEL2_READER_H)
	set(PHOTO2NUMPY_FLAGS -DSKIP_LEVEL2)
endif (NOT HAVE_LEVEL2_READER_H)

i3_add_library(photo2numpy
	private/photo2numpy/photo2numpy.c
	COMPILE_FLAGS ${PHOTO2NUMPY_FLAGS}
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

if(BUILD_PHOTOSPLINE_FITTER)
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
