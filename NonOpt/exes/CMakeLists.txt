# Copyright (C) 2026 Frank E. Curtis
#
# This code is published under the MIT License.
#
# Builds one executable per .cpp in this directory (e.g. solveProblem,
# runExperiment).

file(GLOB NONOPT_EXE_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(exe_src IN LISTS NONOPT_EXE_SOURCES)
  get_filename_component(exe_name "${exe_src}" NAME_WE)
  add_executable(${exe_name} "${exe_src}")
  # CMAKE_DL_LIBS is "-ldl" on Linux and empty on macOS/Windows, so this is
  # portable (the original Makefiles hard-coded -ldl).
  target_link_libraries(${exe_name} PRIVATE NonOptProblems NonOpt ${CMAKE_DL_LIBS})
endforeach()
