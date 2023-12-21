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
add_custom_target(${PROJECT_NAME}-all)
add_custom_target(${PROJECT_NAME}-examples)
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/examples FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    add_dependencies(${PROJECT_NAME}-examples ${tarname})
    add_dependencies(${PROJECT_NAME}-all ${tarname})
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

# add benchmark files
add_custom_target(${PROJECT_NAME}-benchmarks)
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/benchmarks FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    add_dependencies(${PROJECT_NAME}-benchmarks ${tarname})
    add_dependencies(${PROJECT_NAME}-all ${tarname})
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

# add solver files
add_custom_target(${PROJECT_NAME}-solvers)
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR}/solvers FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}")
    set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
    add_dependencies(${PROJECT_NAME}-solvers ${tarname})
    add_dependencies(${PROJECT_NAME}-all ${tarname})
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})

# # add unittests
# aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
# set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
#   CACHE INTERNAL "gismo list of unittests")

# needed?: set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)


