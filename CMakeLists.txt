cmake_minimum_required (VERSION 2.6)

set(CMAKE_BUILD_TYPE Debug CACHE STRING "This is the Debug mode")
project (proffit)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
set(CMAKE_C_FLAGS "-std=c11 -Wall -Wextra -Wno-unused-parameter")


# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}  $ENV{ROOTSYS}/etc/cmake)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net Minuit)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#include(${ROOT_USE_FILE})

FIND_PATH(CFITSIO_INCLUDE_DIR fitsio.h
  HINTS ENV HEADAS
  PATH_SUFFIXES include cfitsio include/cfitsio
  )

FIND_LIBRARY(CFITSIO_LIBRARY cfitsio
  HINTS ENV HEADAS
  PATH_SUFFIXES lib
  )

list(APPEND CMAKE_PREFIX_PATH $PWD)

FIND_PATH(WCS_INCLUDE_DIR wcslib/wcs.h
  HINTS ENV HEADAS
  PATH_SUFFIXES include wcs include/wcs
  )

FIND_LIBRARY(WCS_LIBRARIES wcs PATH $ENV{HEADAS}/lib)

if (WCS_INCLUDE_DIR AND WCS_LIBRARIES)
    set (WCS_FOUND true)
endif(WCS_INCLUDE_DIR AND WCS_LIBRARIES)

if (WCS_FOUND)
    if (NOT WCS_FIND_QUIETLY)
        message(" -- Found wcs: ${WCS_LIBRARIES}")
    endif(NOT WCS_FIND_QUIETLY)
else (WCS_FOUND)
   if (WCS_FIND_REQUIRED)
      message("Fatal_error: Could not find wcs")
    endif(WCS_FIND_REQUIRED)
endif(WCS_FOUND) 

# Find MultiNest library
FIND_PACKAGE(LAPACK REQUIRED)
SET(_INCLUDES ${LAPACK_INCLUDES})
SET(_LIBS ${LAPACK_LIBRARIES})

IF(MULTINEST_INCLUDES)
  # Already in cache, be silent
  set (MULTINEST_FIND_QUIETLY TRUE)
ENDIF (MULTINEST_INCLUDES)

FIND_PATH(MULTINEST_ROOT_DIR
    NAMES include/multinest.h
    HINTS /usr/local/multinest $ENV{MULTINEST_ROOT} 
    DOC "MULTINEST root directory.")
    
FIND_PATH(_MULTINEST_INCLUDE_DIRS
    NAMES multinest.h
    HINTS ${MULTINEST_ROOT_DIR}/include
    DOC "MULTINEST Include directory")
    
FIND_PATH(MULTINEST_MODULE_DIRS
    NAMES nested.mod
    HINTS ${MULTINEST_ROOT_DIR}/modules
    DOC "MULTINEST Fortran module directory")

# Now find the library:
SET(_MULTINEST_LIB_NAME "multinest")
IF(MULTINEST_USE_MPI)
    SET(_MULTINEST_LIB_NAME "multinest_mpi")
ENDIF(MULTINEST_USE_MPI)

FIND_LIBRARY(_MULTINEST_LIBRARY
    NAMES ${_MULTINEST_LIB_NAME}
    HINTS ${MULTINEST_ROOT_DIR}/lib)

SET(MULTINEST_INCLUDE_DIRS ${_MULTINEST_INCLUDE_DIRS} ${_INCLUDES})
SET(MULTINEST_LIBRARIES ${_MULTINEST_LIBRARY} ${_LIBS})

# handle the QUIETLY and REQUIRED arguments and set MULTINEST_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE (FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MULTINEST DEFAULT_MSG MULTINEST_LIBRARIES MULTINEST_INCLUDE_DIRS)

MARK_AS_ADVANCED(MULTINEST_LIBRARIES MULTINEST_INCLUDE_DIRS)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

include_directories(${CFITSIO_INCLUDE_DIR} ${WCS_INCLUDE_DIR} ${ROOT_INCLUDE_DIRS} ${MULTINEST_INCLUDE_DIRS})


add_executable(proffit proffit.cpp)
target_link_libraries(proffit ${CFITSIO_LIBRARY} ${WCS_LIBRARIES} ${ROOT_LIBRARIES} ${MULTINEST_LIBRARIES})



