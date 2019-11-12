# - Try to find BARS installation
# Once done, this will define
#
#  BARS_FOUND     - system has BARS
#  BARSSYS        - build files for BARS
#  BARS_INCLUDES  - the BARS include directories
#  BARS_LIBRARIES - location of libmars.so
#  Author: Alexander Avrorin - avrorin@inr.ru

include(LibFindMacros)


# Paths
#-----------------------------------

find_library(BARS_LIBRARY
  NAMES libmars.so
  HINTS $ENV{BARSSYS}/lib
  NO_DEFAULT_PATH
)

set(BARSSYS $ENV{BARSSYS} CACHE PATH "BARS build directory")
set(BARS_INCLUDE_DIR ${BARSSYS}/include CACHE PATH "BARS include directory")

# Prevent BARS directories from showing up as compilation options
mark_as_advanced(BARS_INCLUDE_DIR BARS_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BARS DEFAULT_MSG BARS_LIBRARY BARS_INCLUDE_DIR BARSSYS)

set (BARS_LIBRARIES ${BARS_LIBRARY})
set (BARS_INCLUDES  ${BARS_INCLUDE_DIR})
