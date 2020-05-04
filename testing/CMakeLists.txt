cmake_minimum_required(VERSION 3.6)

message(STATUS "ACCESSING TEST DIRECTORY: ${PROJECT_SOURCE_DIR}")

# adding gtest subdir with its own cmake file
if(WIN32)
  set(BUILD_SHARED_LIBS ON)
  add_definitions(-DGTEST_HAS_TR1_TUPLE=0)
endif()
add_subdirectory(lib/gtest-1.10.0)

# include gtest library. gtest_SOURCE_DIR is availble from subdir addition
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}/include/gtest ${gtest_SOURCE_DIR})

add_custom_target(unittests)

# Set location of test scripts, test data, and run directory
file(COPY test_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(TEST_SCRIPTS ${CMAKE_CURRENT_SOURCE_DIR}/test_scripts)
set(TEST_DATA ${CMAKE_CURRENT_BINARY_DIR}/test_data)

# Set environment variables for testing.
if(WIN32)
  # Note: Do not indent. The indentation whitespace is added to the variable!
  set(TEST_ENV_PATH_OPTIONS
      "PATH=${CMAKE_CURRENT_BINARY_DIR};\
${CMAKE_BINARY_DIR};\
${CMAKE_BINARY_DIR}/contrib/madlib-1.3.0;\
${CMAKE_BINARY_DIR}/contrib/madlib-1.3.0/Contrib/ANN;\
${gtest_BINARY_DIR};\
${VTK_LIB_DIR};\
${GMSH_LIB_PATH};\
${NETGEN_LIBRARY_DIR};${NETGEN_LIBRARY_DIR}/../bin;\
${ZLIB_LIB_PATH};${ZLIB_LIB_PATH}/../bin;\
${METIS_LIB_PATH};\
${CGNS_LIB_PATH};${CGNS_LIB_PATH}/../bin;\
${HDF5_LIB_PATH};${HDF5_LIB_PATH}/../bin;\
${EXO_LIB_PATH};${EXO_LIB_PATH}/../bin;\
${netCDF_LIB_DIR};${netCDF_LIB_DIR}/../bin;\
${SIMMETRIX_LIB_DIR};\
$ENV{PATH}")
  set(TEST_ENV_LD_OPTIONS
      "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}")
else()
  set(TEST_ENV_PATH_OPTIONS "PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PATH}")
  # Note: Do not indent. The indentation whitespace is added to the variable!
  set(TEST_ENV_LD_OPTIONS
      "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:\
${CMAKE_BINARY_DIR}:\
${VTK_LIB_DIR}:\
${NETGEN_LIBRARY_DIR}:\
${OpenCASCADE_LIBRARY_DIR}:\
$ENV{LD_LIBRARY_PATH}")
endif()

message(STATUS "Test Env Options: ${TEST_ENV_PATH_OPTIONS} ${TEST_ENV_LD_OPTIONS}")

#-------------------------------------------------------------------------------
# NEM_add_test_executable(<test_exe_name>)
#
# Adds a test executable named run<test_exe_name>Test. Requires:
#  ${TEST_SCRIPTS}/test<test_exe_name>.C  = The executable to turn into test
# Generates:
#  run<test_exe_name>Test                 = name of executable
macro(CLX_add_test_executable TEST_EXE_NAME)
  add_executable(run${TEST_EXE_NAME}Test ${TEST_SCRIPTS}/test${TEST_EXE_NAME}.C)
  target_link_libraries(run${TEST_EXE_NAME}Test gtest gtest_main clcx clcxcsc)
  add_dependencies(unittests run${TEST_EXE_NAME}Test)
endmacro()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# NEM_add_test(<test_name> <test_exe_name> <test_pwd> [arguments...])
#
# Adds a test to Nemosys named <test_name>Test using run<test_exe_name>Test
# executable. Requires:
#  run<test_exe_name>Test         = name of executable
#                                   (use NEM_add_test_executable())
#  ${TEST_DATA}/<test_pwd>        = Directory to use as working directory
#                                   to run the test
#  [arguments...]                 = Arguments passed to executable
# Generates:
#  <test_name>Test                = name of CTest test
macro(CLX_add_test TEST_NAME TEST_EXE_NAME TEST_PWD)
  add_test(NAME ${TEST_NAME}Test
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      run${TEST_EXE_NAME}Test ${ARGN}
      WORKING_DIRECTORY ${TEST_DATA}/${TEST_PWD}
  )
endmacro()
#-------------------------------------------------------------------------------

# standard tests
if (ENABLE_CSC)
    CLX_add_test_executable(Beam)
    CLX_add_test(beam Beam beam beam beam-ref.dat)
endif()

