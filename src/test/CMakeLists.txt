project(simcoflowtest Fortran)
set (CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(
    ${Simcoflow_INCLUDE_DIR}
    )
set(HYPRE_INCLUDE_DIR /usr/include/hypre/)
set(${PROJECT_NAME}_INCLUDE_DIR)
set(HYPRE_LIBRARIES /usr/lib/x86_64-linux-gnu/libHYPRE_core.so /usr/lib/x86_64-linux-gnu/libHYPRE.so /usr/lib/x86_64-linux-gnu/libHYPRE_FEI.so)
set (LINK_LIBS ${LINK_LIBS})
macro(my_add_test test)
    add_test(NAME ${test} WORKING_DIRECTORY ${PROJECT_BINARY_DIR} COMMAND
        ${test})
endmacro()
# aux_source_directory(../ SRC_LIST)
list(APPEND SRC_LIST2 
    ../Mesh_Mod.f90
    ../PrecisionVar_Mod.f90)
list(APPEND POINTLIST
    ${SRC_LIST2}
    pointtest.f90)
add_executable(testpoint ${POINTLIST})
target_link_libraries(testpoint ${LINK_LIBS})
my_add_test(testpoint)
list(APPEND CELLLIST
    ${SRC_LIST2}
    celltest.f90)
add_executable(testcell ${CELLLIST})
target_link_libraries(testcell ${LINK_LIBS})
my_add_test(testcell)
list(APPEND PARTICLELIST
    ${SRC_LIST2}
    ../BoundaryFunction.f90
    ../BoundaryFunction2.f90
    ../BoundaryInterface2.f90
    ../Particles.f90
    ../StateVariables.f90
    ../Constants.f90
    ../Clsvof_Mod.f90
    particletest.f90)
add_executable(testparticle ${PARTICLELIST})
target_link_libraries(testparticle ${LINK_LIBS})
my_add_test(testparticle)
list(APPEND PARTICLELIST2
    ${SRC_LIST2}
    ../BoundaryFunction.f90
    ../BoundaryFunction2.f90
    ../BoundaryInterface2.f90
    ../Particles.f90
    ../StateVariables.f90
    ../Constants.f90
    ../Clsvof_Mod.f90
    particletest2.f90)	
add_executable(testparticle2 ${PARTICLELIST2})
target_link_libraries(testparticle2 ${LINK_LIBS})
my_add_test(testparticle2)
list(APPEND CLSVOFTESTLIST
    ${SRC_LIST2}
    ../BoundaryFunction.f90
    ../BoundaryFunction2.f90
    ../BoundaryInterface2.f90
    ../StateVariables.f90
    ../Constants.f90
    ../Cutcell_Mod.f90
    ../Clsvof_Mod.f90
    clsvoftest.f90)	
add_executable(clsvoftest ${CLSVOFTESTLIST})
target_link_libraries(clsvoftest ${LINK_LIBS})
my_add_test(clsvoftest)
