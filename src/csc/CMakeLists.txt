set (CLCXMODSRC
    ${CMAKE_SOURCE_DIR}/src/csc/clcxInterface.C
    ${CMAKE_SOURCE_DIR}/src/csc/clcxCSC.C)
set (CLCXDRVSRC
    ${CMAKE_SOURCE_DIR}/src/csc/clcxDRV.C)

###############
# Targets
###############

# csc module library
add_library(clcxcsc ${CLCXMODSRC})
target_link_libraries(clcxcsc
  PUBLIC
  clcx
  ${MPI_LIBRARIES}
  # ${LAPACK_LIBRARIES}
  # ${BLAS_LIBRARIES}
  ${SPOOLES_LIB}
  ${ARPACK_LIB}
  ${FLANN_LIB}
  ${EXODUS_LIBRARIES}
  ${NetCDF_LIBRARIES}
  Threads::Threads
  IMPACT::SITCOM
  Boost::serialization
  Boost::iostreams
)

target_include_directories(clcxcsc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/clcx>
        Boost::headers
)

target_include_directories(clcxcsc
    SYSTEM PUBLIC
    ${SPOOLES_INC_PATH})

install(TARGETS clcxcsc
    EXPORT clcxcsc
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# csc module driver
add_executable(clcx_drv ${CLCXDRVSRC})
target_include_directories(clcx_drv
    PUBLIC
    IMPACT::SITCOM
    IMPACT::SimOUT)

target_link_libraries(clcx_drv
  ${MPI_LIBRARIES}
  IMPACT::SITCOM
  IMPACT::SimOUT)

install(TARGETS clcx_drv
    RUNTIME DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

