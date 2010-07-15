#
# $Id: CMakeLists.txt 64301 2010-05-25 16:49:40Z jvansanten $
#

i3_project(photospline)

i3_add_library(photospline
	lib/*.c
	USE_TOOLS cfitsio
)
