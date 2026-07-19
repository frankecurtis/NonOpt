# Copyright (C) 2026 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Builds the core NonOpt static library (libNonOpt).

# CONFIGURE_DEPENDS makes CMake re-glob when files are added/removed, so you do
# not have to re-run cmake by hand after adding a new source file.
file(GLOB NONOPT_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB NONOPT_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

add_library(NonOpt STATIC ${NONOPT_SOURCES})
add_library(NonOpt::NonOpt ALIAS NonOpt)

target_include_directories(NonOpt PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/NonOpt>)

# LAPACK first, then BLAS (LAPACK depends on BLAS symbols).
target_link_libraries(NonOpt PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

set_target_properties(NonOpt PROPERTIES
  OUTPUT_NAME NonOpt
  VERSION ${PROJECT_VERSION})

# ---- Install ----------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS NonOpt
        EXPORT NonOptTargets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(FILES ${NONOPT_HEADERS}
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/NonOpt")
