cmake_minimum_required(VERSION 4.3)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

option(BUILD_DEPS OFF)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED
  ${SKBUILD_SABI_COMPONENT}
  )

include(UseCython)
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if (NOT "${SKBUILD_SABI_VERSION}" STREQUAL "")
  if ("${SKBUILD_SABI_VERSION}" VERSION_LESS "3.12")
    message(FATAL_ERROR "ABI3 requires python>=3.12")
  endif()
  set(USE_SABI USE_SABI ${SKBUILD_SABI_VERSION})
endif()

cython_transpile(src/cell_tree2d/cell_tree2d.pyx LANGUAGE CXX OUTPUT_VARIABLE cell_tree2d_c)

python_add_library(cell_tree2d MODULE
                   "${cell_tree2d_c}"
                   src/cell_tree2d_cpp/cell_tree2d_c.h
                   src/cell_tree2d_cpp/cell_tree2d_c.cpp
                   WITH_SOABI
                   ${USE_SABI}
                   )

target_include_directories(cell_tree2d PUBLIC src/cell_tree2d_cpp)

target_link_libraries(cell_tree2d PUBLIC Python::NumPy)

install(TARGETS cell_tree2d DESTINATION cell_tree2d)

