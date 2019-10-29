Program celltest
   USE Mesh, ONLY : cell
   IMPLICIT NONE

   TYPE(cell):: c1
!   c1 = cell(3,4) ! This should not really work, as Isize and Jsize are not set
      call exit(1) ! Not implemented, so should fail
END PROGRAM celltest

