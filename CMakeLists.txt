cmake_minimum_required(VERSION 3.11)
cmake_policy(VERSION 3.11)
include(CPack)
project(simcoflow
    VERSION 0.1
    LANGUAGES Fortran
    )
add_subdirectory(src)
