# Assumed now that all f90 files in same folder, and that there is a folder "under" containing the makefile (e.g. call it rundir)
# In this folder, type make, and files will be compiled and an executable made, here called SIMCOFlow (run code, type: ./SIMCOFlow)
# To delete all mod files etc and executable, type
# make clean 

# since you run on your own pc, you have to change options line below, to where hypre and mpi headerfiles are located
# see p.72 in hypre manual
# https://computation.llnl.gov/sites/default/files/public/hypre-2.11.2_usr_manual.pdf

compiler = mpi
# only mpif77 tested

ifeq ($(compiler),mpi)
	f90 = mpif77
        options = -g -w -I/usr/local/include -I/usr/lib/openmpi/include 
        link_options = -fopenmp /usr/local/lib/libHYPRE.a     
else
	f90 = gfortran
        options = -g -w -I/usr/local/include -I/usr/lib/openmpi/include 
        link_options = -fopenmp /usr/local/lib/libHYPRE.a
endif



objects = Clsvof_Mod.o\
          ComputePUV_Mod.o\
          Cutcell_Mod.o\
          Main.o\
          Matrix_Mod.o\
          Mesh_Mod.o\
          MPI_Mod.o\
          Particles.o\
          PrecisionVar_Mod.o\
          PredictorUV_mod.o\
          PrintResult.o\
          ProjectionP_Mod.o\
          Solver_Mod.o\
          StateVariables.o\
          VTK_mod.o\
          VTR_mod.o 

SIMCOFlow : $(objects)
	$(f90) -o SIMCOFlow $(objects) $(link_options)
#
# -c : lager objekt fil, lenker ikke



Clsvof_Mod.o:../Clsvof_Mod.f90\
				PrecisionVar_Mod.o\
				Cutcell_Mod.o\
				Mesh_Mod.o\
				StateVariables.o\
				Matrix_Mod.o
	$(f90) -c ../Clsvof_Mod.f90 $(options)
#
Cutcell_Mod.o:../Cutcell_Mod.f90 \
				PrecisionVar_Mod.o\
				StateVariables.o\
				Mesh_Mod.o
	$(f90) -c ../Cutcell_Mod.f90 $(options)
#

ComputePUV_Mod.o:../ComputePUV_Mod.f90\
			PrecisionVar_Mod.o\
			Mesh_Mod.o\
			Cutcell_Mod.o\
			Clsvof_Mod.o\
			StateVariables.o\
			PredictorUV_mod.o\
			Particles.o\
			ProjectionP_Mod.o
	$(f90) -c ../ComputePUV_Mod.f90 $(options)
#
VTK_mod.o:../VTK_mod.f90\
			PrecisionVar_Mod.o
	 $(f90) -c ../VTK_mod.f90 $(options)
#
PrintResult.o:../PrintResult.f90 \
			PrecisionVar_Mod.o\
			Mesh_Mod.o\
			Cutcell_Mod.o\
			StateVariables.o\
			VTK_mod.o\
			VTR_mod.o\
			Particles.o\
			Clsvof_Mod.o
	$(f90) -c ../PrintResult.f90 $(options)
#
Main.o:../Main.f90 \
			PrecisionVar_Mod.o\
			Mesh_Mod.o\
			Cutcell_Mod.o\
			StateVariables.o\
			Clsvof_Mod.o\
			StateVariables.o\
			Solver_Mod.o\
			MPI_Mod.o\
			Particles.o\
			PrintResult.o
	$(f90) -c ../Main.f90 $(options)
# 
Mesh_Mod.o:../Mesh_Mod.f90 \
			PrecisionVar_Mod.o\
			MPI_Mod.o
	 $(f90) -c ../Mesh_Mod.f90 $(options)
#
MPI_Mod.o:../MPI_Mod.f90 
	 $(f90) -c ../MPI_Mod.f90 $(options)
# MPI_mod har inclide mpi.h and HYPREf.h....glytt -I hit
PrecisionVar_Mod.o:../PrecisionVar_Mod.f90 
	$(f90) -c ../PrecisionVar_Mod.f90 $(options)
Matrix_Mod.o:../Matrix_Mod.f90 \
			PrecisionVar_Mod.o
	$(f90) -c ../Matrix_Mod.f90 $(options)
#
Particles.o:../Particles.f90 \
						PrecisionVar_Mod.o\
						MPI_Mod.o\
						Mesh_Mod.o\
						StateVariables.o\
						Clsvof_Mod.o
	$(f90) -c ../Particles.f90 $(options)
#
PredictorUV_mod.o:../PredictorUV_mod.f90 \
			PrecisionVar_Mod.o\
			Mesh_Mod.o\
			Clsvof_Mod.o\
			Cutcell_Mod.o\
			StateVariables.o\
			MPI_Mod.o\
			Matrix_Mod.o\
			Particles.o\
			PrintResult.o
	$(f90) -c ../PredictorUV_mod.f90 $(options)
#
StateVariables.o:../StateVariables.f90 \
			PrecisionVar_Mod.o\
			Matrix_Mod.o\
			Mesh_Mod.o
	$(f90) -c ../StateVariables.f90 $(options)
#
VTR_mod.o:../VTR_mod.f90
	$(f90) -c ../VTR_mod.f90 $(options)

#
ProjectionP_Mod.o:../ProjectionP_Mod.f90 \
			PrecisionVar_Mod.o\
			Mesh_Mod.o\
			Cutcell_Mod.o\
			StateVariables.o\
			PredictorUV_mod.o\
			MPI_Mod.o
	$(f90) -c ../ProjectionP_Mod.f90 $(options)
#
Solver_Mod.o : ../Solver_Mod.f90 \
			PrecisionVar_Mod.o\
			Mesh_Mod.o\
			Cutcell_Mod.o\
			Clsvof_Mod.o\
			StateVariables.o\
			ComputePUV_Mod.o\
			MPI_Mod.o\
			Particles.o\
			PrintResult.o

	$(f90) -c ../Solver_Mod.f90 $(options)
#
clean :
	rm *.mod *.o $(objects) SIMCOFlow -f

