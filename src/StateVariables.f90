Module StateVariables
    USE PrecisionVar
    USE Constants, ONLY : g, epsi, roa, row
    USE Mesh, ONLY : TsimcoMesh, Grid, Cell, getMeshSizes
    IMPLICIT NONE
    PRIVATE
!   for run again from
    LOGICAL :: RunAgain
    LOGICAL :: ICorProb
    INTEGER(it8b) :: IttRun
!   For Wave Braking on vertical wall
    TYPE :: TWave
     REAL(KIND=dp) :: t0,cw0,Amp0,Depthw,Lamdaw,twp,HChannel,              &
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
    PUBLIC :: setSolverVariables, getSolverVariables
    contains

    TYPE(TVariables) function constructVar() result( this )
      !
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize, ight, jght
      call getMeshSizes(ibeg, jbeg, Isizee=Isize, Jsizee=Jsize, ighte=ight, jghte=jght)
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
      INTEGER(kind=it4b):: i,j, ibeg, jbeg, Isize, Jsize
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
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
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
   !   call Boundary_Condition_Var(PGrid,PCell,this,0.d0)
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
    subroutine Boundary_Condition_Var(PGrid,PCell,Vari,wave,Time)
      TYPE(Grid),INTENT(IN):: PGrid
      TYPE(Cell),INTENT(IN):: PCell
      TYPE(TVariables),INTENT(INOUT):: Vari
      TYPE(Twave), INTENT(in) :: wave
      REAL(KIND=dp),INTENT(IN):: Time
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize, ight, jght
      INTEGER(kind=it4b):: i,j,temp
      REAL(KIND=dp),PARAMETER:: Twall = 300.d0
      REAL(KIND=dp):: ywu,ywv,xwu,xwv,ywuout,ywvout,xwuout,xwvout,            &
                      Hwin,Hwout,Hain,Haout
      call getMeshSizes(ibeg, jbeg, Isize, Jsize, ight, jght)
      temp=0
      do j=1,Jsize
        if(PCell%vof(1,j)>epsi) Hwin=PGrid%y(1,j)-                             &
                                           (0.5d0-PCell%vof(1,j))*PGrid%dy(1,j)
        if(PCell%vof(Isize,j)>epsi.and.PCell%vof(Isize,j)<1.d0-epsi) then
          Hwout=PGrid%y(Isize,j)-(0.5d0-PCell%vof(Isize,j))*PGrid%dy(Isize,j)
          temp=temp+1
        end if
      end do
      xwu=PGrid%x(1,1)-PGrid%dx(1,1)/2.d0
      xwv=PGrid%x(1,1)-PGrid%dx(1,1)
      ywu=wave%Amp0*dsin(wave%kw*(xwu-wave%cw0*Time))
      ywv=wave%Amp0*dsin(wave%kw*(xwv-wave%cw0*Time))

      xwuout=PGrid%x(Isize,1)+PGrid%dx(Isize,1)/2.d0
      xwvout=PGrid%x(Isize,1)+PGrid%dx(Isize,1)
      ywuout=wave%Amp0*dsin(wave%kw*(xwuout-wave%cw0*Time))
      ywvout=wave%Amp0*dsin(wave%kw*(xwvout-wave%cw0*Time))

      if(temp==2) Hwout=wave%Depthw
      Hain=wave%HChannel-Hwin
      Haout=wave%HChannel-Hwout
      do j = jbeg,Jsize+jbeg-1
      ! Left inlet water
        Vari%p(ibeg-ight,j)=Vari%p(ibeg,j) !Vari%Pint/(Vari%Pref)
        Vari%t(ibeg-ight,j)=Vari%t(ibeg,j)
        if(PGrid%y(1,j)-Hwin<ywu) then
          Vari%u(ibeg-ight,j)=Vari%UwInlet-wave%Amp0*wave%kw*(Vari%UwInlet-wave%cw0)*                   &
                    dsin(wave%kw*(xwu-wave%cw0*Time))*dcosh(wave%kw*PGrid%y(1,j))/dsinh(wave%kw*Hwin)
        else
          Vari%u(ibeg-ight,j)=Vari%UgInlet+wave%Amp0*wave%kw*(Vari%UgInlet-wave%cw0)*                   &
          dsin(wave%kw*(xwu-wave%cw0*Time))*dcosh(wave%kw*(PGrid%y(1,j)-wave%HChannel))/dsinh(wave%kw*Hain)
        end if

        if(PGrid%y(1,j)-Hwin<ywv) then
          Vari%v(ibeg-ight,j)=wave%Amp0*wave%kw*(Vari%UwInlet-wave%cw0)*dcos(wave%kw*(xwv-wave%cw0*Time))*   &
                  dsinh(wave%kw*(PGrid%y(1,j)+0.5d0*PGrid%dy(1,j)))/dsinh(wave%kw*Hwin)
        else
          Vari%v(ibeg-ight,j)=-wave%Amp0*wave%kw*(Vari%UgInlet-wave%cw0)*dcos(wave%kw*(xwv-wave%cw0*Time))*  &
              dsinh(wave%kw*(PGrid%y(1,j)+0.5d0*PGrid%dy(1,j)-wave%HChannel))/dsinh(wave%kw*Hain)
        end if

      ! Right outlet water
      !  if(PGrid%y(Isize,j)-Hwout<ywuout) then
        if(PCell%vof(isize,j)>0.5d0) then
          Vari%u(Isize,j)=Vari%UwInlet!-Amp0*kw*(Vari%UwInlet-cw0)*                       &
            ! dsin(kw*(xwuout-cw0*Time))*dcosh(kw*PGrid%y(Isize,j))/dsinh(kw*Hwout)
        else
          Vari%u(Isize,j)=Vari%UgInlet  !+Amp0*kw*(Vari%UgInlet-cw0)*                       &
            ! dsin(kw*(xwuout-cw0*Time))*dcosh(kw*(PGrid%y(Isize,j)-HChannel))/ &
            !                                                     dsinh(kw*Haout)
        end if

        if(PGrid%y(Isize,j)-Hwout<ywvout) then
          Vari%v(Isize+1,j)=wave%Amp0*wave%kw*(Vari%UwInlet-wave%cw0)*dcos(wave%kw*(xwvout-wave%cw0*Time))*  &
              dsinh(wave%kw*(PGrid%y(Isize,j)+0.5d0*PGrid%dy(Isize,j)))/dsinh(wave%kw*Hwout)
        else
          Vari%v(Isize+1,j)=-wave%Amp0*wave%kw*(Vari%UgInlet-wave%cw0)*dcos(wave%kw*(xwvout-wave%cw0*Time))* &
              dsinh(wave%kw*(PGrid%y(Isize,j)+0.5d0*PGrid%dy(Isize,j)-wave%HChannel))/   &
                                                                  dsinh(wave%kw*Haout)
        end if
        Vari%p(Isize+ibeg+ight-1,j)=Vari%p(Isize+ibeg-1,j)
        Vari%u(Isize+ibeg+ight-1,j)=Vari%u(Isize+ibeg-1,j) !Vari%Uint/Vari%Uref
        Vari%v(Isize+1,j)=Vari%v(Isize+ibeg-1,j)
      !  Vari%v(Isize+1,j)=0.d0
      end do
      do i = ibeg,Isize+ibeg-1
      ! Bottom boundary
        Vari%p(i,jbeg-jght)=Vari%p(i,jbeg)
        Vari%t(i,jbeg-jght)=Vari%t(i,jbeg)
        Vari%u(i,jbeg-jght)=Vari%u(i,jbeg)
        xwu=PGrid%x(i,1)+PGrid%dx(i,1)/2.d0
        Vari%u(i,jbeg-jght)=-wave%Amp0*wave%kw*(Vari%UwInlet-wave%cw0)*dsin(wave%kw*(xwu-wave%cw0*Time))*    &
                        dcosh(wave%kw*(PGrid%y(i,1)-PGrid%dy(i,1)))/dsinh(wave%kw*wave%Depthw)
        xwv=PGrid%x(i,1)
        Vari%v(i,jbeg-jght)=wave%Amp0*wave%kw*(Vari%UwInlet-wave%cw0)*dcos(wave%kw*(xwv-wave%cw0*Time))*    &
                  dsinh(wave%kw*(PGrid%y(i,1)-0.5d0*PGrid%dy(i,1)))/dsinh(wave%kw*wave%Depthw) !0.d0 !Vari%Vint/Vari%Uref
      ! Open air
        Vari%p(i,Jsize+jbeg+jght-1)=0.d0
        Vari%t(i,jbeg+Jsize+jght-1)=Vari%t(i,jbeg+Jsize-1)
        Vari%u(i,Jsize+jbeg+jght-1)=Vari%u(i,Jsize+jbeg-1)
     !   Vari%v(i,Jsize+jbeg-1) = Vari%v(i,Jsize+jbeg-2)
        Vari%v(i,Jsize+jght) = Vari%v(i,Jsize)
      end do
    end subroutine Boundary_Condition_Var

    !The naming scheme is not correct here, FIXME
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


