######################################################################
## CMakeLists.txt --- gsStructuralAnalysis/benchmarks
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris, Hugo Verhelst
## Copyright (C) 2023
######################################################################

# add benchmark files
aux_cpp_directory(${CMAKE_CURRENT_SOURCE_DIR} FILES)
foreach(file ${FILES})
    add_gismo_executable(${file})
    get_filename_component(tarname ${file} NAME_WE) # name without extension
    set_property(TEST ${tarname} PROPERTY LABELS "${PROJECT_NAME}-benchmarks")
    if(GISMO_BUILD_EXAMPLES)
      set_target_properties(${tarname} PROPERTIES FOLDER "${PROJECT_NAME}")
      add_dependencies(${PROJECT_NAME}-benchmarks ${tarname})
      add_dependencies(${PROJECT_NAME}-all ${tarname})
    else(GISMO_BUILD_EXAMPLES)
      set_target_properties(${tarname} PROPERTIES
        FOLDER "${PROJECT_NAME}"
        EXCLUDE_FROM_ALL TRUE)
    endif(GISMO_BUILD_EXAMPLES)
    # install the example executables (optionally)
    install(TARGETS ${tarname} DESTINATION "${BIN_INSTALL_DIR}" COMPONENT exe OPTIONAL)
endforeach(file ${FILES})