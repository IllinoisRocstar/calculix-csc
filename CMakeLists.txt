cmake_minimum_required(VERSION 3.9)

#########################
# project global settings
#########################
set(CLCX_MAJOR_VERSION 2)
set(CLCX_MINOR_VERSION 15)
set(CLCX_PATCH_VERSION 2)
set(CLCX_EXTRA_VERSION "" CACHE STRING "with CSC module")
set(CLCX_VERSION "${CLCX_MAJOR_VERSION}.${CLCX_MINOR_VERSION}")
set(CLCX_VERSION "${CLCX_VERSION}.${CLCX_PATCH_VERSION}${CLCX_EXTRA_VERSION}")
set(CLCX_SHORT_LICENSE "GPL")
set(PROJECT_NAME "CLCX version ${CLCX_VERSION}")
set(CLCX_AUTHOR
    "Masoud Safdari (build system and CSC extensions)"
    "Akash Patel (homogenization extensiosn)"
    "Guido Dhondt (original CalculiX developer)"
    "Klaus Wittig (original CalculiX developer)"
    "Peter A. Gustafson (CalculiX patches)")
set(CLCX_AUTHOR_EMAIL "msafdari@illinoisrocstar.com")
project(${PROJECT_NAME} Fortran CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build.")
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    message(STATUS "System is " ${CMAKE_SYSTEM_NAME})
    set(CMAKE_MACOSX_RPATH 1)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release"
      CACHE STRING "Build type (Release, Debug, RelWithDebugInfo, MinSizeRel)")
endif()

# package module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/)

# set options
option(ENABLE_TESTING "Enable testing" ON)
option(ENABLE_PAR "Builds parallel Epic solver" ON)
option(ENABLE_CSC "Builds solver driver and csc modules" ON)
option(ENABLE_CFD_SUPPORT "Enables usage of offline CFD data" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries (and not static ones)" ON)

# default definitions
add_definitions("-DUSE_MT")
add_definitions("-DARCH=Linux")
#add_definitions("-DMATRIXSTORAGE")
add_definitions("-DEXODUSII")
add_definitions("-DSPOOLES")
add_definitions("-DARPACK")

# setting RPATH
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
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

# According to old_install_manifest, LAPACK and BLAS are only needed as transitive dependencies of TAUCS
# find_package(LAPACK REQUIRED)
# find_package(BLAS REQUIRED)
find_package(Threads REQUIRED)
find_package(Exodus REQUIRED)
find_package(NetCDF REQUIRED)

include_directories(${EXODUS_INCLUDE_DIR})
include_directories(${NetCDF_INCLUDE_DIR})

if (ENABLE_CSC)
  find_package(IMPACT REQUIRED)
  add_definitions("-DHAVE_CSC")
  find_package(Boost 1.48 REQUIRED COMPONENTS serialization iostreams)
endif()

find_library(SPOOLES_LIB spooles spooles.a HINTS /usr/lib /usr/local/lib ENV SPOOLES_DIR)
if (SPOOLES_LIB MATCHES SPOOLES_LIB-NOTFOUND)
    message(FATAL_ERROR "SPOOLES was not found!")
endif()
if (SPOOLES_LIB)
    get_filename_component(SPOOLES_LIB_PATH ${SPOOLES_LIB} DIRECTORY)
    find_path(SPOOLES_INC_PATH SPOOLES.h HINTS ${SPOOLES_LIB_PATH} PATH_SUFFIXES spooles)
    include_directories(${SPOOLES_INC_PATH})
endif()

find_library(ARPACK_LIB arpack arpack_x64 arpack_MAC HINTS /usr/lib /usr/local/lib ENV ARPACK_DIR)
if (ARPACK_LIB MATCHES ARPACK_LIB-NOTFOUND)
    message(FATAL_ERROR "ARPACK was not found!")
endif()

if(ENABLE_CFD_SUPPORT)
    find_library(FLANN_LIB flann PATH_SUFFIXES lib)
    find_path(FLANN_INC "flann.h" PATH_SUFFIXES src include flann)
    if (FLANN_LIB MATCHES FLANN_LIB-NOTFOUND)
        message(FATAL_ERROR "The flann library was not found!")
    endif()
    include_directories(${FLANN_INC})
    add_definitions("-DHAVE_CFD")
endif()

################
# setting files
################
include(GNUInstallDirs)
# install header files
file(GLOB H_FILES include/*.H)
file(GLOB h_FILES include/*.h)
file(GLOB hpp_FILES include/*.hpp)
set(CLCX_INC_FILES ${H_FILES} ${h_FILES} ${hpp_FILES})
install(FILES ${CLCX_INC_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/clcx)

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

#####################
# perform testing
#####################
if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(testing)
endif(ENABLE_TESTING)

############
# packaging
############
set(CPACK_PACKAGE_NAME "Calculix")
# set(CPACK_PACKAGE_VENDOR "")
set(CPACK_PACKAGE_VERSION "${CLCX_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${CLCX_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${CLCX_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${CLCX_PATCH_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Calculix with IMPACT interface")
set(CPACK_PACKAGE_CONTACT "Illinois Rocstar LLC <rocstar@illinoisrocstar.com>")

# Debian
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEBIAN_PACKAGE_RELEASE 1)
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
# set(CPACK_DEBIAN_PACKAGE_DEBUG ON)  # For testing only
set(CPACK_DEBIAN_PACKAGE_CONTROL_STRICT_PERMISSION TRUE)  # CMake 3.4
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)  # CMake 3.6
set(CPACK_DEBIAN_ARCHIVE_TYPE "gnutar")
# RPM
set(CPACK_RPM_FILE_NAME RPM-DEFAULT)
set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_RPM_PACKAGE_RELEASE_DIST ON)
set(CPACK_RPM_PACKAGE_AUTOREQPROV ON)
include(CPack)
