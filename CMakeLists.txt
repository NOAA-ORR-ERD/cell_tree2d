cmake_minimum_required(VERSION 3.24...4.0)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

option(BUILD_DEPS OFF)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

include(UseCython)
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

cython_transpile(src/cell_tree2d/cell_tree2d.pyx LANGUAGE CXX OUTPUT_VARIABLE cell_tree2d_c)

python_add_library(cell_tree2d MODULE
                   "${cell_tree2d_c}"
                   src/cell_tree2d_cpp/cell_tree2d_c.h
                   src/cell_tree2d_cpp/cell_tree2d_c.cpp
                   WITH_SOABI)

target_include_directories(cell_tree2d PUBLIC src/cell_tree2d_cpp)

target_link_libraries(cell_tree2d PUBLIC Python::NumPy)

install(TARGETS cell_tree2d DESTINATION cell_tree2d)

