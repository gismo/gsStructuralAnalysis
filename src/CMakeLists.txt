######################################################################
## CMakeLists.txt --- gsStructuralAnalysis/src
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris, Hugo Verhelst
## Copyright (C) 2023
######################################################################

include_directories(${GISMO_INCLUDE_DIRS})

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})

if (GISMO_BUILD_PCH)
  add_precompiled_header(gsPrecompiledHeader misc/gsPrecompiledHeader.h)
endif()

  ## Add module
  add_library(gsStructuralAnalysis_pybind OBJECT
    gsStructuralAnalysis_.cpp
    )
  set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:gsStructuralAnalysis_pybind>
    CACHE INTERNAL "G+Smo modules" )

FOREACH(subdir ${SUBDIRS})

  get_filename_component(GM_NAME ${subdir} NAME)
  #message("module: ${GM_NAME}")

  ## Collect files
  aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_H  )
  aux_cpp_noins_directory  (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_CPP)
  aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_HPP)

  if(GISMO_BUILD_LIB)
    aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR}/${subdir} ${subdir}_INS)
  endif()

  ## Add module
  add_library(${GM_NAME} OBJECT
    ${${subdir}_H}
    ${${subdir}_HPP}
    ${${subdir}_CPP} # static/non templated part
    ${${subdir}_INS}
    )

  set_target_properties(${GM_NAME} PROPERTIES
    COMPILE_DEFINITIONS gismo_EXPORTS
    POSITION_INDEPENDENT_CODE ON
    LINKER_LANGUAGE CXX
    FOLDER "G+Smo modules"
    )

  if(TARGET gmp)
    add_dependencies(${GM_NAME} gmp)
  endif()
  if(TARGET mpfr)
    add_dependencies(${GM_NAME} mpfr)
  endif()

  add_dependencies(${GM_NAME} ${PROJECT_NAME})

  ## Add extra instances (template instantiations)
  # NOTE: The static/non templated part of the library is
  # compiled above only for GISMO_COEFF_TYPE
  math(EXPR ii "1")
  foreach(ins ${GISMO_EXTRA_INSTANCE})
    add_library(${GM_NAME}_${ii} OBJECT
      ${${subdir}_H}
      ${${subdir}_HPP}
      ${${subdir}_INS}
      )
    set_target_properties(${GM_NAME}_${ii} PROPERTIES
      COMPILE_DEFINITIONS "gismo_EXPORTS;real_t=${ins}"
      POSITION_INDEPENDENT_CODE ON
      LINKER_LANGUAGE CXX
      FOLDER "G+Smo modules"
      )
    set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${GM_NAME}_${ii}>
      CACHE INTERNAL "G+Smo modules" )
    math(EXPR ii "${ii}+1")
  endforeach()

  set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${GM_NAME}>
    CACHE INTERNAL "G+Smo modules" )

ENDFOREACH()
