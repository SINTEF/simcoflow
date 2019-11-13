Module StateVariables
    USE PrecisionVar
    USE BoundaryFunction2, ONLY : BCBase2
    USE Constants, ONLY : g, epsi, roa, row
    USE Mesh, ONLY : TsimcoMesh, Grid, Cell, Isize, Jsize, ibeg, &
         &           jbeg, ight, jght
    USE BoundaryFunction, ONLY : BCBase, BCUW, BCUE, BCUS, BCUN, BCVW,         &
                                 BCVE, BCVS, BCVN, BCPW, BCPE, BCPS, BCPN
    IMPLICIT NONE
    PRIVATE
!   for run again from
    LOGICAL :: RunAgain
    LOGICAL :: ICorProb
    INTEGER(it8b) :: IttRun
!   For Wave Braking on vertical wall
    TYPE :: TWave
      REAL(KIND=dp) :: t0,cw0,Amp0,Depthw,Lamdaw,twp,HChannel,                 &
                      LChannel,kw,omew,HDomain
    END TYPE TWave
    TYPE :: TVariables
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: u,v,p,t,Gpu,Gpv,ures,vres,pres,mres
      REAL(KIND=dp):: Uint,Vint,Pint,Tint,Uref,Roref,Pref,Tref
      REAL(dp) :: UwInlet, UgInlet, Rey, Fr, nuref
    CONTAINS
      PROCEDURE, PASS(thiS), PUBLIC :: Initialize
    END TYPE TVariables
    interface Boundary_Condition_Var
      module procedure Boundary_Condition_Var
    end interface Boundary_Condition_Var
    interface TVariables
       module procedure constructVar
    end interface
    interface TWave
       module procedure constructWave
    end interface
    Interface getSolverVariables
       Module Procedure getSolverVariables
    End Interface getSolverVariables
    Interface setSolverVariables
       Module Procedure setSolverVariables
    End Interface setSolverVariables
    PUBLIC :: TWave, TVariables
    PUBLIC:: Boundary_Condition_Var
    PUBLIC:: Boundary_Condition_Var2
    PUBLIC :: setSolverVariables, getSolverVariables
    contains

    TYPE(TVariables) function constructVar() result( this )
      !
      allocate(this%u(ibeg-ight:Isize-ibeg+ight+1,jbeg-jght:Jsize-jbeg+jght+1))
      allocate(this%v(ibeg-ight:Isize-ibeg+ight+1,jbeg-jght:Jsize-jbeg+jght+1))
      allocate(this%p(ibeg-ight:Isize-ibeg+ight+1,jbeg-jght:Jsize-jbeg+jght+1))
      allocate(this%Gpu(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
      allocate(this%Gpv(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
      allocate(this%ures(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
      allocate(this%vres(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
      allocate(this%pres(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
      allocate(this%mres(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
      allocate(this%t(ibeg-ight:Isize-ibeg+ight+1,jbeg-jght:Jsize-jbeg+jght+1))
    end function constructVar

    subroutine Initialize(this, wave, simcomesh, PCell,PGrid,Uint,Vint,Pint,&
         &                Tint,Uref,Tref,Roref,Lref, nuref, ha, UwInlet, UgInlet )
      class(TVariables),INTENT(INOUT) :: this
      TYPE(TWave), INTENT(in)         :: wave
      TYPE(TsimcoMesh), INTENT(in) :: simcomesh
      REAL(KIND=dp),INTENT(IN):: Uint,Vint,Pint,Tint,Uref,Tref,Roref,Lref, nuref,ha
      REAL(dp),  INTENT(in) :: UwInlet, UgInlet
      TYPE(Cell),INTENT(IN):: PCell
      TYPE(Grid),INTENT(IN):: PGrid
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp):: Hwt0,xu0,Amp,cw,tw,delta,xwu,xwv,ywu,ywv
      this%Uint = Uint
      this%Vint = Vint
      this%Pint = Pint
      this%Tint = Tint
      this%Uref = Uref
      this%Roref = Roref
      this%Pref = Roref*Uref**2.d0
      this%Tref = Tref
      this%nuref = nuref
      this%UwInlet = UwInlet
      this%UgInlet = UgInlet
      do i = ibeg,Isize+ibeg-1
        do j = jbeg,Jsize+jbeg-1
        ! for sinusoidal wave
          this%u(i,j)=0.d0
          this%v(i,j)=0.d0
          Hwt0=wave%Depthw
          if(PGrid%y(i,j)>=Hwt0) then
            this%p(i,j)=roa/roref*g*(PGrid%y(1,Jsize)+PGrid%dy(1,Jsize)/2.d0-  &
                                                                  PGrid%y(i,j))
          else
            this%p(i,j)=roa/Roref*g*(PGrid%y(1,Jsize)+PGrid%dy(1,Jsize)/2.d0-  &
              Hwt0)+row/Roref*g*(Hwt0-PGrid%y(i,j))
          endif
          this%Gpu(i,j)=0.d0
          this%Gpv(i,j)=0.d0
          this%t(i,j)=Tint/Tref
        end do
      end do
      do i=1,Isize
        do j=1,Jsize
          xwu=PGrid%x(i,j)+PGrid%dx(i,j)/2.d0
          xwv=PGrid%x(i,j)
          ywu=wave%Amp0*dsin(wave%kw*xwu)
          ywv=wave%Amp0*dsin(wave%kw*xwv)
          if(PGrid%y(i,j)-wave%Depthw<ywu) then
            this%u(i,j)=UwInlet-wave%Amp0*wave%kw*(UwInlet-wave%cw0)*dsin(wave%kw*xwu)*            &
                        dcosh(wave%kw*PGrid%y(i,j))/dsinh(wave%kw*wave%Depthw)
          else
            this%u(i,j)=UgInlet+wave%Amp0*wave%kw*(UgInlet-wave%cw0)*dsin(wave%kw*xwu)*            &
                               dcosh(wave%kw*(PGrid%y(i,j)-wave%HChannel))/dsinh(wave%kw*Ha)
          end if
          if(PGrid%y(i,j)-wave%Depthw<ywv) then
            this%v(i,j)=wave%Amp0*wave%kw*(UwInlet-wave%cw0)*dcos(wave%kw*xwv)*                    &
               dsinh(wave%kw*(PGrid%y(i,j)+0.5d0*PGrid%dy(i,j)))/dsinh(wave%kw*wave%Depthw)
          else
            this%v(i,j)=-wave%Amp0*wave%kw*(UgInlet-wave%cw0)*dcos(wave%kw*xwv)*                   &
               dsinh(wave%kw*(PGrid%y(i,j)+0.5d0*PGrid%dy(i,j)-wave%HChannel))/dsinh(wave%kw*Ha)
          end if
        end do
      end do
      this%Rey=Uref*Lref/nuref
      this%Fr=Uref/dsqrt(g*Lref)
      print*,"Reynolds number:",this%Rey
      print*,"Froude number:",this%Fr
    end subroutine Initialize

    ! This is kind of artificial
    TYPE(TWave) function constructWave(t0, cw0, Amp0, Depthw, Lamdaw, twp, HChannel, &
         &                             LChannel, kw, omew, HDomain ) RESULT ( this )
      REAL(dp), INTENT(in) :: t0
      REAL(dp), INTENT(in) :: cw0
      REAL(dp), INTENT(in) :: Amp0
      REAL(dp), INTENT(in) :: Depthw
      REAL(dp), INTENT(in) :: Lamdaw
      REAL(dp), INTENT(in) :: twp
      REAL(dp), INTENT(in) :: HChannel
      REAL(dp), INTENT(in) :: LChannel
      REAL(dp), INTENT(in) :: kw
      REAL(dp), INTENT(in) :: omew
      REAL(dp), INTENT(in) :: HDomain

      this % t0 = t0
      this % cw0 = cw0
      this % Amp0 = Amp0
      this % Depthw = Depthw
      this % Lamdaw = Lamdaw
      this % twp = twp
      this % HChannel = HChannel
      this % LChannel = LChannel
      this % kw = kw
      this % omew = omew
      this % HDomain = HDomain
    end function constructWave

    !*******************************************************
    !
    !
    !      ________Slip Wall_______
    ! Inlet
    ! Air
    !      ________________________
    ! Inlet........................Outlet
    ! Water........................Water
    !      ........................
    !      ........................
    !      __________Slip Wall_____
    !
    !*******************************************************
    subroutine Boundary_Condition_Var(PGrid, PCell, Vari, BCp, BCu, BCv, Time)
      TYPE(Grid), INTENT(IN)  	      :: PGrid
      TYPE(Cell), INTENT(IN)  	      :: PCell
      TYPE(TVariables), INTENT(INOUT) :: Vari
      TYPE(BCBase), INTENT(INOUT)     :: BCu, BCv, BCp  	
      REAL(KIND=dp),INTENT(IN)        :: Time
      INTEGER(kind=it4b)              :: i,j

      ! Compute the boundary value for u-velocity.
      ! Western boundary 
      call BCUW(BCu, PGrid%x(1,:)-PGrid%dx(1,:)/2.d0, PGrid%y(1,:), 	       &
                PGrid%dx(1,:), PGrid%dy(1,:), Vari%p(1,:), Vari%u(1,:),        &
                Vari%v(1,:), PCell%vof(1,:), PCell%phi(1,:), Time)
      Vari%u(ibeg-ight,:) = BCu%VarW(:)
      ! Eastern boundary
      call BCUE(BCu, PGrid%x(Isize,:)+PGrid%dx(Isize,:)/2.d0, PGrid%y(Isize,:),&
                PGrid%dx(Isize,:), PGrid%dy(Isize,:), Vari%p(Isize,:),         &
                Vari%u(Isize,:), Vari%v(Isize,:), PCell%vof(Isize,:),          &
                PCell%phi(Isize,:), Time)            
      Vari%u(Isize+ight,:) = BCu%VarE(:)  
      ! Southern boundary
      call BCUS(BCu, PGrid%x(:,1)+PGrid%dx(:,1)/2.d0, 			       &
  	        PGrid%y(:,1)-PGrid%dy(:,1)/2.d0, PGrid%dx(:,1),                &
                PGrid%dy(:,1), Vari%p(:,1), Vari%u(:,1),           	       &
                Vari%v(:,1), PCell%vof(:,1), PCell%phi(:,1), Time)
      if(BCu%flag(3) == 0) then
        Vari%u(:,jbeg-jght) = BCu%VarS(:)-Vari%u(:,1)    
      else
        Vari%u(:,jbeg-jght) = 2.d0*BCu%VarS(:)-Vari%u(:,1)  
      end if
      ! Northern boundary
      call BCUN(BCu, PGrid%x(:,Jsize)+PGrid%dx(:,Jsize)/2.d0, 		       &
  	        PGrid%y(:,Jsize)+PGrid%dy(:,Jsize)/2.d0, PGrid%dx(:,Jsize),    &
                PGrid%dy(:,Jsize), Vari%p(:,Jsize), Vari%u(:,Jsize),           &
                Vari%v(:,Jsize), PCell%vof(:,Jsize), PCell%phi(:,Jsize), Time)
      if(BCu%flag(4) == 0) then
        Vari%u(:,Jsize+jght) = BCu%VarS(:)-Vari%u(:,Jsize)    
      else
        Vari%u(:,Jsize+jght) = 2.d0*BCu%VarS(:)-Vari%u(:,Jsize)  
      end if
      
      ! Compute the boundary value for v-velocity
      ! Western boundary
      call BCVW(BCv, PGrid%x(1,:)-PGrid%dx(1,:)/2.d0,                          &
                PGrid%y(1,:)+PGrid%dy(1,:)/2.d0, PGrid%dx(1,:), 	       &
                PGrid%dy(1,:), Vari%p(1,:), Vari%u(1,:),                       &
                Vari%v(1,:), PCell%vof(1,:), PCell%phi(1,:), Time)  
      if(BCv%flag(1) == 0) then
        Vari%v(ibeg-ight,:) = BCv%VarW(:)-Vari%v(ibeg,:)
      else
        Vari%v(ibeg-ight,:) = 2.d0*BCv%VarW(:)-Vari%v(ibeg,:)
      end if   
      ! Eastern boundary
      call BCVE(BCv, PGrid%x(Isize,:)+PGrid%dx(Isize,:)/2.d0,                  &
                PGrid%y(Isize,:)+PGrid%dy(Isize,:)/2.d0, PGrid%dx(Isize,:),    &
                PGrid%dy(Isize,:), Vari%p(Isize,:), Vari%u(Isize,:),           &
                Vari%v(Isize,:), PCell%vof(Isize,:), PCell%phi(Isize,:), Time)  
      if(BCv%flag(2) == 0) then
        Vari%v(Isize+ight,:) = BCv%VarE(:)-Vari%v(Isize,:)
      else
        Vari%v(Isize+ight,:) = 2.d0*BCv%VarE(:)-Vari%v(Isize,:)
      end if
      ! Southern boundary
      call BCVS(BCv, PGrid%x(:,1), PGrid%y(:,1)-PGrid%dy(:,1)/2.d0, 	       &
                PGrid%dx(:,1), PGrid%dy(:,1), Vari%p(:,1), Vari%u(:,1),        &
                Vari%v(:,1), PCell%vof(:,1), PCell%phi(:,1), Time)
      Vari%v(:,jbeg-jght) = BCv%VarS(:)
      ! Northern boundary 
      call BCVN(BCv, PGrid%x(:,Jsize), PGrid%y(:,Jsize)-PGrid%dy(:,Jsize)/2.d0,&
                PGrid%dx(:,Jsize), PGrid%dy(:,Jsize), Vari%p(:,Jsize),         &
                Vari%u(:,Jsize), Vari%v(:,Jsize), PCell%vof(:,Jsize),  	       &
                PCell%phi(:,Jsize), Time) 
      Vari%v(:,Jsize+jght) = BCv%VarN(:)
      
      ! Compute the boundary value for pressure
      ! Western boundary  
      call BCpW(BCp, PGrid%x(1,:)+PGrid%dx(1,:)/2.d0, PGrid%y(1,:),            &
                PGrid%dx(1,:), PGrid%dy(1,:), Vari%p(1,:), Vari%u(1,:),        &
                Vari%v(1,:), PCell%vof(1,:), PCell%phi(1,:), Time)
      Vari%p(ibeg-ight,:) = BCp%VarW(:)
      ! Eastern boundary 
      call BCpE(BCp, PGrid%x(Isize,:)+PGrid%dx(Isize,:)/2.d0, PGrid%y(Isize,:),&
                PGrid%dx(Isize,:), PGrid%dy(Isize,:), Vari%p(Isize,:),         &
                Vari%u(Isize,:), Vari%v(Isize,:), PCell%vof(Isize,:),           &
                PCell%phi(Isize,:), Time)
      Vari%p(Isize+ight,:) = BCp%VarE(:)
      ! Southern boundary
      call BCpS(BCp, PGrid%x(:,1), PGrid%y(:,1)-PGrid%dy(:,1)/2.d0,            &
                PGrid%dx(:,1), PGrid%dy(:,1), Vari%p(:,1), Vari%u(:,1),        &
                Vari%v(:,1), PCell%vof(:,1), PCell%phi(:,1), Time)
      Vari%p(:,jbeg-jght) = BCp%VarS(:)
      ! Northern boundary 
      call BCpN(BCp, PGrid%x(:,Jsize), PGrid%y(:,Jsize)-PGrid%dy(:,Jsize)/2.d0,&
                PGrid%dx(:,Jsize), PGrid%dy(:,Jsize), Vari%p(:,Jsize),         &
                Vari%u(:,Jsize), Vari%v(:,Jsize), PCell%vof(:,Jsize),           &
                PCell%phi(:,Jsize), Time)
      Vari%p(:,Jsize+jght) = BCp%VarN(:)      
    end subroutine Boundary_Condition_Var

    !*******************************************************
    !
    !
    !      ________Slip Wall_______
    ! Inlet
    ! Air
    !      ________________________
    ! Inlet........................Outlet
    ! Water........................Water
    !      ........................
    !      ........................
    !      __________Slip Wall_____
    !
    !*******************************************************
    subroutine Boundary_Condition_Var2(PGrid, PCell, Vari, BCp, BCu, BCv, Time)
      TYPE(Grid), INTENT(IN)  	      :: PGrid
      TYPE(Cell), INTENT(IN)  	      :: PCell
      TYPE(TVariables), INTENT(INOUT) :: Vari
      TYPE(BCBase2), INTENT(INOUT)     :: BCu, BCv, BCp  	
      REAL(KIND=dp),INTENT(IN)        :: Time
      INTEGER(kind=it4b)              :: i,j

      ! Compute the boundary value for u-velocity.
      ! Western boundary 
      call BCu%west(PGrid%x(1,:)-PGrid%dx(1,:)/2.d0, PGrid%y(1,:), 	       &
                PGrid%dx(1,:), PGrid%dy(1,:), Vari%p(1,:), Vari%u(1,:),        &
                Vari%v(1,:), PCell%vof(1,:), PCell%phi(1,:), Time)
      Vari%u(ibeg-ight,:) = BCu%VarW(:)
      ! Eastern boundary
      call BCU%east(PGrid%x(Isize,:)+PGrid%dx(Isize,:)/2.d0, PGrid%y(Isize,:),&
                PGrid%dx(Isize,:), PGrid%dy(Isize,:), Vari%p(Isize,:),         &
                Vari%u(Isize,:), Vari%v(Isize,:), PCell%vof(Isize,:),          &
                PCell%phi(Isize,:), Time)            
      Vari%u(Isize+ight,:) = BCu%VarE(:)  
      ! Southern boundary
      call BCU%south(PGrid%x(:,1)+PGrid%dx(:,1)/2.d0, 			       &
  	        PGrid%y(:,1)-PGrid%dy(:,1)/2.d0, PGrid%dx(:,1),                &
                PGrid%dy(:,1), Vari%p(:,1), Vari%u(:,1),           	       &
                Vari%v(:,1), PCell%vof(:,1), PCell%phi(:,1), Time)
      if(BCu%flag(3) == 0) then
        Vari%u(:,jbeg-jght) = BCu%VarS(:)-Vari%u(:,1)    
      else
        Vari%u(:,jbeg-jght) = 2.d0*BCu%VarS(:)-Vari%u(:,1)  
      end if
      ! Northern boundary
      call BCU%north(PGrid%x(:,Jsize)+PGrid%dx(:,Jsize)/2.d0, 		       &
  	        PGrid%y(:,Jsize)+PGrid%dy(:,Jsize)/2.d0, PGrid%dx(:,Jsize),    &
                PGrid%dy(:,Jsize), Vari%p(:,Jsize), Vari%u(:,Jsize),           &
                Vari%v(:,Jsize), PCell%vof(:,Jsize), PCell%phi(:,Jsize), Time)
      if(BCu%flag(4) == 0) then
        Vari%u(:,Jsize+jght) = BCu%VarS(:)-Vari%u(:,Jsize)    
      else
        Vari%u(:,Jsize+jght) = 2.d0*BCu%VarS(:)-Vari%u(:,Jsize)  
      end if
      
      ! Compute the boundary value for v-velocity
      ! Western boundary
      call BCV%west(PGrid%x(1,:)-PGrid%dx(1,:)/2.d0,                          &
                PGrid%y(1,:)+PGrid%dy(1,:)/2.d0, PGrid%dx(1,:), 	       &
                PGrid%dy(1,:), Vari%p(1,:), Vari%u(1,:),                       &
                Vari%v(1,:), PCell%vof(1,:), PCell%phi(1,:), Time)  
      if(BCv%flag(1) == 0) then
        Vari%v(ibeg-ight,:) = BCv%VarW(:)-Vari%v(ibeg,:)
      else
        Vari%v(ibeg-ight,:) = 2.d0*BCv%VarW(:)-Vari%v(ibeg,:)
      end if   
      ! Eastern boundary
      call BCV%east(PGrid%x(Isize,:)+PGrid%dx(Isize,:)/2.d0,                  &
                PGrid%y(Isize,:)+PGrid%dy(Isize,:)/2.d0, PGrid%dx(Isize,:),    &
                PGrid%dy(Isize,:), Vari%p(Isize,:), Vari%u(Isize,:),           &
                Vari%v(Isize,:), PCell%vof(Isize,:), PCell%phi(Isize,:), Time)  
      if(BCv%flag(2) == 0) then
        Vari%v(Isize+ight,:) = BCv%VarE(:)-Vari%v(Isize,:)
      else
        Vari%v(Isize+ight,:) = 2.d0*BCv%VarE(:)-Vari%v(Isize,:)
      end if
      ! Southern boundary
      call BCV%south(PGrid%x(:,1), PGrid%y(:,1)-PGrid%dy(:,1)/2.d0, 	       &
                PGrid%dx(:,1), PGrid%dy(:,1), Vari%p(:,1), Vari%u(:,1),        &
                Vari%v(:,1), PCell%vof(:,1), PCell%phi(:,1), Time)
      Vari%v(:,jbeg-jght) = BCv%VarS(:)
      ! Northern boundary 
      call BCV%north(PGrid%x(:,Jsize), PGrid%y(:,Jsize)-PGrid%dy(:,Jsize)/2.d0,&
                PGrid%dx(:,Jsize), PGrid%dy(:,Jsize), Vari%p(:,Jsize),         &
                Vari%u(:,Jsize), Vari%v(:,Jsize), PCell%vof(:,Jsize),  	       &
                PCell%phi(:,Jsize), Time) 
      Vari%v(:,Jsize+jght) = BCv%VarN(:)
      
      ! Compute the boundary value for pressure
      ! Western boundary  
      call BCp%west(PGrid%x(1,:)+PGrid%dx(1,:)/2.d0, PGrid%y(1,:),            &
                PGrid%dx(1,:), PGrid%dy(1,:), Vari%p(1,:), Vari%u(1,:),        &
                Vari%v(1,:), PCell%vof(1,:), PCell%phi(1,:), Time)
      Vari%p(ibeg-ight,:) = BCp%VarW(:)
      ! Eastern boundary 
      call BCp%east(PGrid%x(Isize,:)+PGrid%dx(Isize,:)/2.d0, PGrid%y(Isize,:),&
                PGrid%dx(Isize,:), PGrid%dy(Isize,:), Vari%p(Isize,:),         &
                Vari%u(Isize,:), Vari%v(Isize,:), PCell%vof(Isize,:),           &
                PCell%phi(Isize,:), Time)
      Vari%p(Isize+ight,:) = BCp%VarE(:)
      ! Southern boundary
      call BCp%south(PGrid%x(:,1), PGrid%y(:,1)-PGrid%dy(:,1)/2.d0,            &
                PGrid%dx(:,1), PGrid%dy(:,1), Vari%p(:,1), Vari%u(:,1),        &
                Vari%v(:,1), PCell%vof(:,1), PCell%phi(:,1), Time)
      Vari%p(:,jbeg-jght) = BCp%VarS(:)
      ! Northern boundary 
      call BCp%north(PGrid%x(:,Jsize), PGrid%y(:,Jsize)-PGrid%dy(:,Jsize)/2.d0,&
                PGrid%dx(:,Jsize), PGrid%dy(:,Jsize), Vari%p(:,Jsize),         &
                Vari%u(:,Jsize), Vari%v(:,Jsize), PCell%vof(:,Jsize),           &
                PCell%phi(:,Jsize), Time)
      Vari%p(:,Jsize+jght) = BCp%VarN(:)      
    end subroutine Boundary_Condition_Var2

    ! The naming scheme is not correct here, FIXME
    SUBROUTINE setSolverVariables(RunAgain_, ICorProb_, IttRun_)
      !
      LOGICAL,       INTENT(in) :: RunAgain_, ICorProb_
      INTEGER(it8b), INTENT(in) :: IttRun_
      !
      RunAgain = RunAgain_
      ICorProb = ICorProb_
      IttRun   = IttRun_
    END SUBROUTINE setSolverVariables

    SUBROUTINE getSolverVariables(RunAgain_, ICorProb_, IttRun_)
      !
      LOGICAL, OPTIONAL, INTENT(out) :: RunAgain_, ICorProb_
      INTEGER(it8b), OPTIONAL, INTENT(out) :: IttRun_
      !
      IF(PRESENT(RunAgain_)) RunAgain_ = RunAgain
      IF(PRESENT(ICorProb_)) ICorProb_ = ICorProb
      IF(PRESENT(IttRun_)) IttRun_     = IttRun
    END SUBROUTINE getSolverVariables

End Module StateVariables


