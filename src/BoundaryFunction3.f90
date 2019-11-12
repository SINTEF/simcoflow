Module BoundaryFunction3
   use BoundaryInterface
  USE PrecisionVar
  USE Constants
   IMPLICIT NONE
   PRIVATE
  PUBLIC :: BCLvsN
  CONTAINS
  SUBROUTINE BCLvsN(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for level set function 
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize
   
    ArrSize = sizeof(this%VarN)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarN, lvsin, -dyin, this%flag(4), this%const(4), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCLvsN

  SUBROUTINE TypicalBC(Arr, Varin, dxy, flag, const, ArrSize)
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: Arr
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)    		    :: Varin, dxy
    INTEGER(KIND=it4b), INTENT(IN)    	                    :: flag, ArrSize
    REAL(KIND=dp), INTENT(IN)				    :: const
    INTEGER(KIND=it4b)	   	  			    :: i
    
    IF(flag == 0) then ! Dirichlet BC
      DO i = 1, ArrSize
        Arr(i) = const
      END DO
    ELSE ! Neumann BC
      DO i = 1, ArrSize
        Arr(i) = Varin(i) - const*dxy(i)/2.d0
      END DO 
    END IF
  END SUBROUTINE TypicalBC   
END MODULE BoundaryFunction3

