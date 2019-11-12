cmake_minimum_required(VERSION 3.4.3 FATAL_ERROR)
get_filename_component (BARS_APP ${CMAKE_CURRENT_SOURCE_DIR} NAME)

project(${BARS_APP})

add_executable(${BARS_APP} 
  src/main.cc 
  src/opts.cc
)

# You need to tell CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
# Set path to CMake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Dependencies
#-----------------------------------

#---Locate the BZip2 package
find_package(BZip2 REQUIRED)
if(NOT BZip2_FOUND)
    message(SEND_ERROR "Failed to find BZip2")
    return()
endif()

#---Locate the ZLIB package
find_package(ZLIB REQUIRED)
if(NOT ZLIB_FOUND)
    message(SEND_ERROR "Failed to find ZLIB")
    return()
endif()

#---Locate GSL
find_package(GSL REQUIRED)

#---Locate CERN ROOT
find_package(ROOT "6" REQUIRED COMPONENTS Gui MathMore TMVA)
if(NOT ROOT_FOUND)
  message(SEND_ERROR "Failed to find ROOT")
  return()
else()
  if(ROOT_VERSION VERSION_LESS "6")
    message(SEND_ERROR "Found ROOT version is ${ROOT_VERSION}, need at least 6")
    return()
  endif()
endif()

#--Locate BARS
find_package(BARS REQUIRED)
if(NOT BARS_FOUND)
  message(SEND_ERROR "Failed to find BARS")
  return()
endif()

# Includes
#-----------------------------------
# Include ROOT
include_directories(${ROOT_INCLUDE_DIRS})
# Include BARS
file(GLOB_RECURSE l  "${BARSSYS}/include" "${BARSSYS}/include/*")
set(list_inc_dirs)
foreach(elem ${l})
  get_filename_component(elem_path ${elem} DIRECTORY)
  list(APPEND list_inc_dirs ${elem_path})
endforeach()
list(REMOVE_DUPLICATES list_inc_dirs)
foreach(elem ${list_inc_dirs})
  include_directories(${elem})
endforeach()


# Compile flags
#-----------------------------------
add_compile_options(-std=c++11)
add_compile_options(-Wall)
add_compile_options(-Winit-self)
add_compile_options(-fPIC)
add_compile_options(-Wpointer-arith)
add_compile_options(-Wcast-align)
add_compile_options(-Woverloaded-virtual)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS})

# Linking
#-----------------------------------
#---The link will fail if there are unresolved symbols
if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
    include_directories("/opt/local/include")
else()
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

target_link_libraries(${BARS_APP} ${BARS_LIBRARIES})
target_link_libraries(${BARS_APP} ${ROOT_LIBRARIES})
target_link_libraries(${BARS_APP} ${BZIP2_LIBRARIES})
target_link_libraries(${BARS_APP} ${ZLIB_LIBRARIES})
target_link_libraries(${BARS_APP} ${GSL_LIBRARIES})

# Installation
#-----------------------------------

# Install binary
install(
  FILES "${CMAKE_BINARY_DIR}/${BARS_APP}"
  DESTINATION  "$ENV{HOME}/.local/bin/"
  PERMISSIONS 
    OWNER_READ 
    OWNER_WRITE 
    OWNER_EXECUTE 
    GROUP_READ 
    GROUP_EXECUTE 
    WORLD_READ 
    WORLD_EXECUTE
)

# Install config
install(
  FILES "${CMAKE_SOURCE_DIR}/src/config.rc"
  DESTINATION  "$ENV{HOME}/.config/bars/default-configs"
  PERMISSIONS 
    OWNER_READ 
    OWNER_WRITE 
    GROUP_READ 
    WORLD_READ 
  RENAME ${BARS_APP}.rc
)
