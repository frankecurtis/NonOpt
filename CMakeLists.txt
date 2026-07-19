# Copyright (C) 2026 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Top-level CMake build for NonOpt.
#
# Place this file at the repository root (next to the existing "NonOpt/" source
# directory, "NonOpt-Manual/", and "docs/"). The per-directory CMakeLists.txt
# files go in NonOpt/src, NonOpt/problems, NonOpt/exes, and NonOpt/tests.
#
# Quick start:
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j
#   ctest --test-dir build --output-on-failure
#   cmake --install build --prefix /your/install/prefix     # optional

cmake_minimum_required(VERSION 3.15)

project(NonOpt
        VERSION 2.1
        DESCRIPTION "NonOpt: Nonconvex, Nonsmooth Optimizer"
        HOMEPAGE_URL "https://github.com/frankecurtis/NonOpt"
        LANGUAGES CXX)

# ---- Default build type -----------------------------------------------------
# A numerical solver should be optimized by default; the original Makefiles
# built with "-g" only (no -O), which can be many times slower.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
               STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# ---- Language settings ------------------------------------------------------
set(CMAKE_CXX_STANDARD 11)            # bump to 17 if/when you modernize
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Allow the static library to be linked into shared objects (Python/AMPL wrappers).
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Put all build artifacts in predictable places.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# ---- Warnings ---------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  add_compile_options(-Wall)
  # Uncomment to tighten up incrementally:
  # add_compile_options(-Wextra -Wpedantic)
endif()

# ---- Dependencies -----------------------------------------------------------
# BLAS and LAPACK are required and are NOT bundled with NonOpt.
#   Ubuntu/Debian : sudo apt-get install libblas-dev liblapack-dev
#   Fedora/RHEL   : sudo dnf install blas-devel lapack-devel
#   macOS         : uses the Accelerate framework automatically,
#                   or: brew install openblas lapack
#
# To force a specific vendor implementation, e.g. OpenBLAS or MKL:
#   cmake -S . -B build -DBLA_VENDOR=OpenBLAS
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# ---- Options ----------------------------------------------------------------
option(NONOPT_BUILD_EXES  "Build the example executables (solveProblem, runExperiment)" ON)
option(NONOPT_BUILD_TESTS "Build the test suite" ON)

# ---- Targets ----------------------------------------------------------------
add_subdirectory(NonOpt/src)
add_subdirectory(NonOpt/problems)

if(NONOPT_BUILD_EXES)
  add_subdirectory(NonOpt/exes)
endif()

if(NONOPT_BUILD_TESTS)
  enable_testing()
  add_subdirectory(NonOpt/tests)
endif()

# ---- Package config export --------------------------------------------------
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(EXPORT NonOptTargets
        NAMESPACE NonOpt::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/NonOpt")

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/NonOptConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion)

# A minimal config file so downstream projects can do find_package(NonOpt).
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/NonOptConfig.cmake.in"
"@PACKAGE_INIT@
include(CMakeFindDependencyMacro)
find_dependency(BLAS)
find_dependency(LAPACK)
include(\"\${CMAKE_CURRENT_LIST_DIR}/NonOptTargets.cmake\")
")
configure_package_config_file(
  "${CMAKE_CURRENT_BINARY_DIR}/NonOptConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/NonOptConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/NonOpt")

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/NonOptConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/NonOptConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/NonOpt")
