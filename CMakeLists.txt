cmake_minimum_required(VERSION 3.6)

#########################
# project global settings
#########################
set(CLCX_MAJOR_VERSION 2)
set(CLCX_MINOR_VERSION 15)
set(CLCX_PATCH_VERSION 0)
set(CLCX_EXTRA_VERSION "" CACHE STRING "with CSC Module")
set(CLCX_VERSION "${CLCX_MAJOR_VERSION}.${CLCX_MINOR_VERSION}")
set(CLCX_VERSION "${CLCX_VERSION}.${CLCX_PATCH_VERSION}${CLCX_EXTRA_VERSION}")
set(CLCX_SHORT_LICENSE "GPL")
set(PROJECT_NAME "CLCX version ${CLCX_VERSION}")

project(${PROJECT_NAME} Fortran CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release"
      CACHE STRING "Build type (Release, Debug, RelWithDebugInfo, MinSizeRel)")
endif()

set(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")
macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro(opt)

# package module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/)

# set options
opt(PAR "Builds parallel Epic solver" ON)
opt(CSC "Builds solver driver and csc modules" ON)

# default definitions
add_definitions("-DUSE_MT")
add_definitions("-DARCH=Linux")
add_definitions("-DMATRIXSTORAGE")
add_definitions("-DSPOOLES")
add_definitions("-DARPACK")
add_definitions("-DEXODUSII")

# compiled objects output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# setting RPATH
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

####################
# finding libraries
####################
if(ENABLE_PAR)
  find_package(MPI REQUIRED)
  # Workaround if MPI is found but MPI_FOUND is not set (windows+msmpi+findmpi)
  set(MPI_FOUND TRUE CACHE BOOL "")
  mark_as_advanced(MPI_FOUND)
  include_directories(${MPI_C_INCLUDE_PATH})
  include_directories(${MPI_CXX_INCLUDE_PATH})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "MPI Libs are " ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_Fortran_LIBRARIES})
  # MPI libraries may be needed for testing
  set(CMAKE_C_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
  set(CMAKE_CXX_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})
  set(CMAKE_Fortran_REQUIRED_LIBRARIES ${MPI_Fortran_LIBRARIES})
  # needed for Calculix
  add_definitions("-DCALCULIX_MPI")
endif(ENABLE_PAR)

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Exodus REQUIRED)
find_package(Threads REQUIRED)
if (ENABLE_CSC)
  find_package(IMPACT REQUIRED)
  add_definitions("-DHAVE_CSC")
endif()

find_library(SPOOLES_LIB spooles HINTS /usr/lib /usr/local/lib)
if (SPOOLES_LIB MATCHES SPOOLES_LIB-NOTFOUND)
    message(FATAL_ERROR "SPOOLES was not found!")
endif()

find_library(ARPACK_LIB arpack HINTS /usr/lib /usr/local/lib)
if (ARPACK_LIB MATCHES ARPACK_LIB-NOTFOUND)
    message(FATAL_ERROR "ARPACK was not found!")
endif()

find_library(EXOII_LIB exoIIv2c HINTS /usr/lib /usr/local/lib)
if (EXOII_LIB MATCHES EXOII_LIB-NOTFOUND)
    message(FATAL_ERROR "EXODUS II was not found!")
endif()

################
# setting files
################
# install header files
file(GLOB H_FILES include/*.H)
file(GLOB h_FILES include/*.h)
file(GLOB hpp_FILES include/*.hpp)
set(CLCX_INC_FILES ${H_FILES} ${h_FILES} ${hpp_FILES})
install(FILES ${CLCX_INC_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

add_subdirectory(src)

#####################
# set soname version
#####################
set_target_properties(clcx PROPERTIES VERSION ${CLCX_VERSION}
  SOVERSION ${CLCX_MAJOR_VERSION})
if (ENABLE_CSC)
  set_target_properties(clcxcsc PROPERTIES VERSION ${CLCX_VERSION}
    SOVERSION ${CLCX_MAJOR_VERSION})
endif()


############
# packaging
############
set(CPACK_PACKAGE_NAME "Calculix")
set(CPACK_PACKAGE_VENDOR "")
set(CPACK_PACKAGE_VERSION "${CLCX_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${CLCX_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${CLCX_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${CLCX_PATCH_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Calculix with IMPACT interface")
set(CPACK_PACKAGE_HOMEPAGE_URL "")
set(CPACK_PACKAGE_CONTACT "msafdari@illinoisrocstar.com")

if(CMAKE_VERSION VERSION_LESS 3.10)
  message(WARNING "CMake 3.10 or above is required to build a DEB package.")
else()
  set(CPACK_GENERATOR DEB)
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "${CPACK_PACKAGE_HOMEPAGE_URL}")
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  set(CPACK_DEBIAN_PACKAGE_DEBUG ON)  # For testing only
  set(CPACK_DEBIAN_PACKAGE_CONTROL_STRICT_PERMISSION TRUE)  # CMake 3.4
  set(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")  # CMake 3.6
  set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)  # CMake 3.6
endif()
include(CPack)

