cmake_minimum_required(VERSION 3.11)
cmake_policy(VERSION 3.11)
include(CPack)
project(simcoflow
    VERSION 0.1
    LANGUAGES Fortran
    )
option(WITH_TESTS "Compile tests" ON)
# Tests
# Run them with "ctest -V" or "make check" to get verbose output.
if (WITH_TESTS)
  enable_testing()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose -F)
endif (WITH_TESTS)
if(WITH_TESTS)
#  add_subdirectory(tests/simcoflow)
endif(WITH_TESTS)
add_subdirectory(src)
