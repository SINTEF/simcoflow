MODULE BoundaryInterface
  USE PrecisionVar

  IMPLICIT NONE
  PRIVATE

  TYPE :: BCBase2
    INTEGER(KIND=it4b)                     :: flag(4)
    !
    REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: Const
    REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: VarW,VarE,VarS,VarN 
    !
    PROCEDURE(SetDN),       PASS(this), PUBLIC, POINTER :: SetDN=>NULL()
    PROCEDURE(SetConstant), PASS(this), PUBLIC, POINTER :: SetConstant=>NULL()
    PROCEDURE(BCinterface), PASS(this), PUBLIC, POINTER :: north=>NULL()
    PROCEDURE(BCinterface), PASS(this), PUBLIC, POINTER :: south=>NULL()
    PROCEDURE(BCinterface), PASS(this), PUBLIC, POINTER :: east=>NULL()
    PROCEDURE(BCinterface), PASS(this), PUBLIC, POINTER :: west=>NULL()
  END TYPE BCBase2
  ABSTRACT INTERFACE
     SUBROUTINE BCinterface(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
        IMPORT :: dp, BCBase2
        Class(BCBase2), INTENT(INOUT)            :: this
        REAL(KIND=dp),  DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin
        REAL(KIND=dp),  DIMENSION(:), INTENT(IN) :: pin, uin, vin
        REAL(KIND=dp),  DIMENSION(:), INTENT(IN) :: vofin, lvsin
        REAL(KIND=dp),                INTENT(IN) :: time 
     END SUBROUTINE BCinterface
!
     SUBROUTINE SetDN(this, W, E, N, S)
       IMPORT it4b
       IMPORT BCBase2
     !! Set the Dirichlet or Neumann boundary condition, 0 : Dirichlet, 1 : Neumann
       CLASS(BCBase2), INTENT(INOUT)      :: this
       INTEGER(KIND=it4b), INTENT(IN)    :: W, E, N, S

     END SUBROUTINE SetDN

     SUBROUTINE SetConstant(this, ConstIn)
        IMPORT :: dp, BCBase2
        CLASS(BCBase2),                           INTENT(INOUT) :: this
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: ConstIn

     END SUBROUTINE SetConstant
  END INTERFACE

  PUBLIC :: BCBase2

END MODULE BoundaryInterface
