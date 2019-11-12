Module BoundaryFunction2
   use BoundaryInterface
  USE PrecisionVar
  USE Constants
   IMPLICIT NONE
   PRIVATE

   INTERFACE BCBase2
      MODULE PROCEDURE construct
   END INTERFACE
   PUBLIC :: BCBase2
  PUBLIC :: BCUW, BCUE, BCUN, BCUS, 						&
            BCVW, BCVE, BCVS, BCVN, 						&
            BCPW, BCPE, BCPS, BCPN, 						&
            BCVofW, BCVofE, BCVofS, BCVofN, 					&
            BCLvsW, BCLvsE, BCLvsS, BCLvsN 		 
   CONTAINS

      type(BCBase2) function construct(Isize, Jsize) RESULT(this)

            INTEGER(KIND=it4b), INTENT(IN) :: Isize, Jsize
            ALLOCATE(this%VarW(Jsize))
            ALLOCATE(this%VarE(Jsize))
            ALLOCATE(this%VarS(Isize))
            ALLOCATE(this%VarN(Isize))
            this%SetConstant => SetConstant
            this%SetDN => SetDN
      end function
          
  SUBROUTINE SetDN(this, W, E, N, S)
  !! Set the Dirichlet or Neumann boundary condition, 0 : Dirichlet, 1 : Neumann    
    CLASS(BCBase2), INTENT(INOUT)      :: this
    INTEGER(KIND=it4b), INTENT(IN)    :: W, E, N, S
    
    this%flag(1) = W
    this%flag(2) = E
    this%flag(3) = S
    this%flag(4) = N
  END SUBROUTINE SetDN
  SUBROUTINE SetConstant(this, ConstIn)
    CLASS(BCBase2), INTENT(INOUT) 		     	 :: this
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: ConstIn	
    
    ALLOCATE(this%Const(sizeof(ConstIn))) 
    this%Const(:) = ConstIn(:) 
  END SUBROUTINE SetConstant
    
  SUBROUTINE BCUW(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for u velocity
    CLASS(BCBase2), INTENT(INOUT)       	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: j, ArrSize
    REAL(KIND=dp)			    :: yw, xw, Hwin, Hain

    ! Western boundary condition 
    ! For simple boundary condition   
    ArrSize = sizeof(this%VarW)
    CALL TypicalBC(this%VarW, uin, dxin, this%flag(1), this%const(1), ArrSize) 
    ! For user defined boundary condition.
    ! In this case, inlet wave boundary condition.    
    DO j = 1, ArrSize
      IF(vofin(j) > epsi) Hwin = yin(j) - (0.5d0 - vofin(j))*dyin(j) ! Compute water level at inlet
    END DO 
    IF(Hwin < epsi) Hwin = this%Const(11)
    Hain = this%Const(1) - Hwin ! Compute the gas height
    xw = xin(1)
    yw = this%Const(5)*dsin(this%Const(6)*(xw-this%Const(7)*time)) 
    DO j = 1, Arrsize
      IF(yin(j)-Hwin < yw) THEN
        this%VarW(j) = this%Const(8) - this%Const(5)*this%Const(6)*            &
                       (this%Const(8) - this%Const(7))*dsin(this%Const(6)*     &
                       (xw - this%Const(7)*Time))*dcosh(this%Const(6)*yin(j))/ &
                       dsinh(this%Const(6)*Hwin)
      ELSE
        this%VarW(j) = this%Const(9) + this%Const(5)*this%Const(6)*            &
                      (this%Const(8) - this%Const(7))*dsin(this%Const(6)*      &
                      (xw - this%Const(7)*Time))*dcosh(this%Const(6)*(yin(j)-  &
                       this%Const(10)))/dsinh(this%Const(6)*Hain)
      END IF
    END DO
  END SUBROUTINE BCUW     
  
  SUBROUTINE BCUE(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for u velocity.
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: ArrSize
    
    ArrSize = sizeof(this%VarE)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarE, uin, -dxin, this%flag(2), this%const(2), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCUE
  
  SUBROUTINE BCUS(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for u velocity.
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: ArrSize
    
    ArrSize = sizeof(this%VarS) 
    ! For simple boundary condition       
    CALL TypicalBC(this%VarS, uin, dyin, this%flag(3), this%const(3), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCUS
  
  SUBROUTINE BCUN(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for u velocity.
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize

    ArrSize = sizeof(this%VarN)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarN, uin, -dyin, this%flag(4), this%const(4), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCUN 
  
  SUBROUTINE BCVW(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for v velocity
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: j, ArrSize
    REAL(KIND=dp)		  	    :: yw, xw, Hwin, Hain
    
    ArrSize = sizeof(this%VarW)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarW, vin, dxin, this%flag(1), this%const(1), ArrSize) 
    ! For user defined boundary condition.
    ! In this case, inlet wave boundary condition.    
    DO j = 1, ArrSize
      IF(vofin(j) > epsi) Hwin = yin(j) - (0.5d0 - vofin(j))*dyin(j) ! Compute water level at inlet
    END DO 
    IF(Hwin < epsi) Hwin = this%Const(11)
    Hain = this%Const(1) - Hwin ! Compute the gas height
    xw = xin(1)
    yw = this%Const(5)*dsin(this%Const(6)*(xw-this%Const(7)*time)) 
    DO j = 1, ArrSize
      IF(yin(j)-Hwin < yw) THEN
        this%VarW(j) = this%Const(5)*this%Const(6)*(this%Const(8) -            &
                       this%Const(7))*dcos(this%Const(6)*                      &
                       (xw - this%Const(7)*Time))*   			       &	
                       dsinh(this%Const(6)*yin(j))/dsinh(this%Const(6)*Hwin)
      ELSE
        this%VarW(j) = this%Const(5)*this%Const(6)*(this%Const(9) -            &
                       this%Const(7))*dcos(this%Const(6)*                      &
                       (xw - this%Const(7)*Time))*                             &
                       dsinh(this%Const(6)*(yin(j) - this%Const(10)))/         &
                       dsinh(this%Const(6)*Hain)
      END IF
    END DO 
  END SUBROUTINE BCVW
  
  SUBROUTINE BCVE(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for v velocity
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: ArrSize
    
    ArrSize = sizeof(this%VarE)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarE, vin, -dxin, this%flag(2), this%const(2), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCVE 
  
  SUBROUTINE BCVS(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for v velocity
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: ArrSize
    
    ArrSize = sizeof(this%VarS)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarS, vin, dyin, this%flag(3), this%const(3), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCVS 
  
  SUBROUTINE BCVN(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for v velocity
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: ArrSize
  
    ArrSize = sizeof(this%VarN)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarN, vin, -dyin, this%flag(4), this%const(4), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCVN
  
  SUBROUTINE BCPW(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for pressure
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)	   	    :: time 
    INTEGER(KIND=it4b)			    :: ArrSize
    
    ArrSize = sizeof(this%VarW)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarW, pin, dxin, this%flag(1), this%const(1), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCPW
  
  SUBROUTINE BCPE(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for pressure
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: ArrSize

    ArrSize = sizeof(this%VarE)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarE, pin, -dxin, this%flag(2), this%const(2), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCPE

  SUBROUTINE BCPS(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for pressure
    CLASS(BCBase2), INTENT(INOUT) 	    :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)		    :: time 
    INTEGER(KIND=it4b)			    :: ArrSize
    
    ArrSize = sizeof(this%VarS)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarS, pin, dyin, this%flag(3), this%const(3), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCPS

  SUBROUTINE BCPN(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for pressure
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize
    
    ArrSize = sizeof(this%VarN)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarN, pin, -dyin, this%flag(4), this%const(4), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCPN
  
  SUBROUTINE BCVofW(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for volume of fluid
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize

    ArrSize = sizeof(this%VarW)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarW, vofin, dxin, this%flag(1), this%const(1), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCVofW
  
  SUBROUTINE BCVofE(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for volume of fluid
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize
    
    ArrSize = sizeof(this%VarE)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarE, vofin, -dxin, this%flag(2), this%const(2), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCVofE

  SUBROUTINE BCVofS(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for volume of fluid
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize

    ArrSize = sizeof(this%VarS)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarS, vofin, dyin, this%flag(3), this%const(3), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCVofS
  
  SUBROUTINE BCVofN(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the northern boundary for volume of fluid
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize
    
    ArrSize = sizeof(this%VarN)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarN, vofin, -dyin, this%flag(4), this%const(4), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCVofN
  
  SUBROUTINE BCLvsW(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the western boundary for level set function
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize
  
    ArrSize = sizeof(this%VarW)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarW, lvsin, dxin, this%flag(1), this%const(1), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCLvsW
  
  SUBROUTINE BCLvsE(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the eastern boundary for level set function 
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize
    
    ArrSize = sizeof(this%VarE) 
    ! For simple boundary condition.   
    CALL TypicalBC(this%VarE, lvsin, -dxin, this%flag(2), this%const(2), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCLvsE

  SUBROUTINE BCLvsS(this, xin, yin, dxin, dyin, pin, uin, vin, vofin, lvsin, time)
  !! Compute the boundary value at the southern boundary for level set function 
    CLASS(BCBase2), INTENT(INOUT) 		   	 :: this
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: xin, yin, dxin, dyin 
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: pin, uin, vin
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: vofin, lvsin
    REAL(KIND=dp), INTENT(IN)				 :: time 
    INTEGER(KIND=it4b)					 :: ArrSize
    
    ArrSize = sizeof(this%VarS)
    ! For simple boundary condition   
    CALL TypicalBC(this%VarS, lvsin, dyin, this%flag(3), this%const(3), ArrSize) 
    ! For user defined boundary condition.   
  END SUBROUTINE BCLvsS
  
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
END MODULE BoundaryFunction2
