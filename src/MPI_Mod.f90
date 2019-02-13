Module MPI
  IMPLICIT NONE
  Include 'mpif.h'
  Include 'HYPREf.h'
  INTEGER:: ierr,num_procs,myid
  Contains
    Subroutine MPI_Initial
      Call MPI_Init(ierr)
      Call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
      Call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    End subroutine MPI_Initial
End module MPI
