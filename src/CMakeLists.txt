project(Simcoflow Fortran)
find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
# Hypre,
# This should be ectracted out to a FindHYPRE.cmake
# The paths here are platform dependent. Please edit
set(HYPRE_INCLUDE_DIR /usr/include/hypre)
set(${PROJECT_NAME}_INCLUDE_DIR)
set(HYPRE_LIBRARIES /usr/lib/libHYPRE.so /usr/lib/libHYPRE_FEI.so /usr/lib/libHYPRE_struct_mv.so /usr/lib/libHYPRE_IJ_mv.so /usr/lib/libHYPRE_parcsr_ls.so)
set(HYPRE_LIBRARIES /usr/lib/x86_64-linux-gnu/libHYPRE_core.so /usr/lib/x86_64-linux-gnu/libHYPRE.so /usr/lib/x86_64-linux-gnu/libHYPRE_FEI.so)
include_directories(${MPI_Fortran_INCLUDE_PATH} ${HYPRE_INCLUDE_DIR} ${Simcoflow_INCLUDE_DIR})
link_directories(${MPI_Fortran_LIBRARIES})
aux_source_directory(. SRC_LIST)
list(APPEND SRC_LIST 
    VTKFortran/src/third_party/BeFoR64/src/lib/befor64.F90
    VTKFortran/src/third_party/BeFoR64/src/lib/befor64_pack_data_m.F90
    VTKFortran/src/third_party/FoXy/src/lib/foxy.f90
    VTKFortran/src/third_party/FoXy/src/lib/foxy_xml_file.f90
    VTKFortran/src/third_party/FoXy/src/lib/foxy_xml_tag.F90
    VTKFortran/src/third_party/PENF/src/lib/penf_b_size.F90
    VTKFortran/src/third_party/PENF/src/lib/penf.F90
    VTKFortran/src/third_party/PENF/src/lib/penf_global_parameters_variables.F90
    VTKFortran/src/third_party/PENF/src/lib/penf_stringify.F90
    VTKFortran/src/third_party/StringiFor/src/lib/stringifor.F90
    VTKFortran/src/third_party/StringiFor/src/lib/stringifor_string_t.F90
    VTKFortran/src/lib/vtk_fortran_dataarray_encoder.f90
    VTKFortran/src/lib/vtk_fortran.f90
    VTKFortran/src/lib/vtk_fortran_parameters.f90
    VTKFortran/src/lib/vtk_fortran_pvtk_file.f90
    VTKFortran/src/lib/vtk_fortran_vtk_file.f90
    VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_abstract.f90
    VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_appended.f90
    VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_ascii_local.f90
    VTKFortran/src/lib/vtk_fortran_vtk_file_xml_writer_binary_local.f90
    VTKFortran/src/lib/vtk_fortran_vtm_file.F90
    Solver_Mod.f90
    PrecisionVar_Mod.f90 
    VTK_mod.f90 
    MPI_Mod.f90 
    Mesh_Mod.f90 
    VTR_mod.f90 
    Cutcell_Mod.f90 
    Clsvof_Mod.f90 
    Particles.f90
    PrintResult.f90 
    PredictorUV_mod.f90
    ProjectionP_Mod.f90 
    ComputePUV_Mod.f90
    BoundaryInterface2.f90
    BoundaryFunction2.f90
    StateVariables.f90 
    Constants.f90
    BoundaryFunction.f90
    BoundaryFunction3.f90
    Main.f90)
set (LINK_LIBS ${LINK_LIBS} 
    ${MPI_Fortran_LIBRARIES}
    ${HYPRE_LIBRARIES}
    )
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})
    
add_subdirectory(test)
