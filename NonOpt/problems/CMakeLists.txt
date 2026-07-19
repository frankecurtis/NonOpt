# Copyright (C) 2026 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Builds the NonOpt test-problem library (libNonOptProblems).

file(GLOB NONOPT_PROBLEM_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(NonOptProblems STATIC ${NONOPT_PROBLEM_SOURCES})
add_library(NonOpt::Problems ALIAS NonOptProblems)

target_include_directories(NonOptProblems PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# Inherits the src include directory and BLAS/LAPACK via the NonOpt target.
target_link_libraries(NonOptProblems PUBLIC NonOpt)

include(GNUInstallDirs)
install(TARGETS NonOptProblems
        EXPORT NonOptTargets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
