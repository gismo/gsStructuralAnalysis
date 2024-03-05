######################################################################
## CMakeLists.txt --- gsStructuralAnalysis
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris, Hugo Verhelst
## Copyright (C) 2023
######################################################################

## gsStructuralAnalysis module
project(gsStructuralAnalysis)

# Apply G+Smo config
include(gsConfig)

set (GISMO_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src"
     CACHE INTERNAL "${PROJECT_NAME} include directories" )

# Add main sources directory
add_subdirectory(src)

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_INS)
  if(${PROJECT_NAME}_INS)
    LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
  endif()
endif()

# add filedata folder
add_definitions(-DSTRAN_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/filedata/")

# add example files
if(GISMO_BUILD_EXAMPLES)
    add_custom_target(${PROJECT_NAME}-all)
    add_custom_target(${PROJECT_NAME}-examples)
    add_custom_target(${PROJECT_NAME}-benchmarks)
    add_custom_target(${PROJECT_NAME}-solvers)
    add_custom_target(${PROJECT_NAME}-tutorials)
    add_subdirectory(examples)
    add_subdirectory(benchmarks)
    add_subdirectory(solvers)
    add_subdirectory(tutorials)
else()
    add_subdirectory(examples EXCLUDE_FROM_ALL)
    add_subdirectory(benchmarks EXCLUDE_FROM_ALL)
    add_subdirectory(solvers EXCLUDE_FROM_ALL)
    add_subdirectory(tutorials EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_EXAMPLES)

# # add unittests
# aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
# set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
#   CACHE INTERNAL "gismo list of unittests")

# needed?: set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)


