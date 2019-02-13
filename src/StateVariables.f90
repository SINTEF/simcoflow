Module StateVariables
    USE PrecisionVar
    USE Mesh
    USE Matrix
    IMPLICIT NONE
    PRIVATE
    INTEGER(kind=it4b),PUBLIC::ight=1,jght=1
    INTEGER(kind=it4b),PUBLIC::ite
    TYPE(Point),PUBLIC:: Start_Point,End_Point
    REAL(KIND=dp),PARAMETER,PUBLIC::pi=4.d0*datan(1.d0),Cp=1.005d3,            &
                                    kT=0.0271d0,kTw=0.0271d0,g = 9.81d0,factor=0.5d0
    REAL(KIND=dp),PUBLIC::Lref,Rey,Fr,wa,Ta,xc,yc,UwInlet,UgInlet,             &
                          xmax,Roref,Hw,Ha,VofInlet,nuref
!   for run again from
    INTEGER(KIND=it8b),PUBLIC::IttRun,IttBegin
    LOGICAL,PUBLIC:: RunAgain
!   For particles from 3D to 2D
    REAL(KIND=dp),PUBLIC::zp,UParInlet,HParInlet,DParInlet,rop,gx,gy
!   For number of particles into system
    INTEGER(kind=it4b),PUBLIC:: NParInlet,IParInlet
!   For Water Entry problem
    REAL(KIND=dp),PUBLIC::y0,y1
!   For Wave Braking on vertical wall
    REAL(KIND=dp),PUBLIC::t0,cw0,Amp0,Depthw,Lamdaw,twp,HChannel,              &
                          LChannel,kw,omew
    REAL(KIND=dp),PUBLIC,PARAMETER::nuw=1.0034d-6,nua=1.506d-5,                &
                                    roa=1.225d0,row=998.2d0 ! At T = 20oC
    REAL(KIND=dp),PUBLIC,PARAMETER::epsi=1.d-3,epsiF=1.d-2,BetaVis=0.5d0
    CHARACTER*70,PUBLIC,PARAMETER:: dir='/home/sontd/code/BoomCase/Data/'
    TYPE,PUBLIC  :: Variables
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: u,v,p,t,Gpu,Gpv,ures,vres,pres,mres
      REAL(KIND=dp):: Uint,Vint,Pint,Tint,Uref,Roref,Pref,Tref
    END TYPE
    PUBLIC:: Initial_Var,Boundary_Condition_Var
    interface Initial_Var
      module procedure Initial_Var
    end interface Initial_Var
    interface Boundary_Condition_Var
      module procedure Boundary_Condition_Var
    end interface Boundary_Condition_Var
    contains
    subroutine Initial_Var(PCell,PGrid,Vari,Uint,Vint,Pint,Tint,Uref,Tref,Roref,Lref)
      REAL(KIND=dp),INTENT(IN):: Uint,Vint,Pint,Tint,Uref,Tref,Roref,Lref
      TYPE(Cell),INTENT(IN):: PCell
      TYPE(Grid),INTENT(IN):: PGrid
      TYPE(Variables),INTENT(INOUT):: Vari
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp):: Hwt0,xu0,Amp,cw,tw,delta,xwu,xwv,ywu,ywv
      Vari%Uint = Uint
      Vari%Vint = Vint
      Vari%Pint = Pint
      Vari%Tint = Tint
      Vari%Uref = Uref
      Vari%Roref = Roref
      Vari%Pref = Roref*Uref**2.d0
      Vari%Tref = Tref
      do i = ibeg,Isize+ibeg-1
        do j = jbeg,Jsize+jbeg-1
        ! for sinusoidal wave
          Vari%u(i,j)=0.d0
          Vari%v(i,j)=0.d0
          Hwt0=Depthw
          if(PGrid%y(i,j)>=Hwt0) then
            Vari%p(i,j)=roa/roref*g*(PGrid%y(1,Jsize)+PGrid%dy(1,Jsize)/2.d0-  &
                                                                  PGrid%y(i,j))
          else
            Vari%p(i,j)=roa/Roref*g*(PGrid%y(1,Jsize)+PGrid%dy(1,Jsize)/2.d0-  &
              Hwt0)+row/Roref*g*(Hwt0-PGrid%y(i,j))
          endif
          Vari%Gpu(i,j)=0.d0
          Vari%Gpv(i,j)=0.d0
          Vari%t(i,j)=Tint/Tref
        end do
      end do
      do i=1,Isize
        do j=1,Jsize
          xwu=PGrid%x(i,j)+PGrid%dx(i,j)/2.d0
          xwv=PGrid%x(i,j)
          ywu=Amp0*dsin(kw*xwu)
          ywv=Amp0*dsin(kw*xwv)
          if(PGrid%y(i,j)-Hw<ywu) then
            Vari%u(i,j)=UwInlet-Amp0*kw*(UwInlet-cw0)*dsin(kw*xwu)*            &
                                          dcosh(kw*PGrid%y(i,j))/dsinh(kw*Hw)
          else
            Vari%u(i,j)=UgInlet+Amp0*kw*(UgInlet-cw0)*dsin(kw*xwu)*            &
                               dcosh(kw*(PGrid%y(i,j)-HChannel))/dsinh(kw*Ha)
          end if
          if(PGrid%y(i,j)-Hw<ywv) then
            Vari%v(i,j)=Amp0*kw*(UwInlet-cw0)*dcos(kw*xwv)*                    &
                        dsinh(kw*(PGrid%y(i,j)+0.5d0*PGrid%dy(i,j)))/dsinh(kw*Hw)
          else
            Vari%v(i,j)=-Amp0*kw*(UgInlet-cw0)*dcos(kw*xwv)*                   &
               dsinh(kw*(PGrid%y(i,j)+0.5d0*PGrid%dy(i,j)-HChannel))/dsinh(kw*Ha)
          end if
        end do
      end do
      Rey=Uref*Lref/nuref
      Fr=Uref/dsqrt(g*Lref)
      print*,"Reynolds number:",Rey
      print*,"Froude number:",Fr
   !   call Boundary_Condition_Var(PGrid,Vari,0.d0)
    end subroutine Initial_Var

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
    subroutine Boundary_Condition_Var(PGrid,Vari,Time)
      TYPE(Grid),INTENT(IN):: PGrid
      TYPE(Variables),INTENT(INOUT):: Vari
      REAL(KIND=dp),INTENT(IN):: Time
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp),PARAMETER:: Twall = 300.d0
      REAL(KIND=dp):: ywu,ywv,xwu,xwv
      xwu=PGrid%x(1,1)-PGrid%dx(1,1)/2.d0
      xwv=PGrid%x(1,1)-PGrid%dx(1,1)
      ywu=Amp0*dsin(kw*(xwu-cw0*Time))
      ywv=Amp0*dsin(kw*(xwv-cw0*Time))

      do j = jbeg,Jsize+jbeg-1
      ! Left inlet water
        Vari%p(ibeg-ight,j)=Vari%p(ibeg,j) !Vari%Pint/(Vari%Pref)
        Vari%t(ibeg-ight,j)=Vari%t(ibeg,j)
        if(PGrid%y(1,j)-Hw<ywu) then
          Vari%u(ibeg-ight,j)=UwInlet-Amp0*kw*(UwInlet-cw0)*                   &
                    dsin(kw*(xwu-cw0*Time))*dcosh(kw*PGrid%y(1,j))/dsinh(kw*Hw)
        else
          Vari%u(ibeg-ight,j)=UgInlet+Amp0*kw*(UgInlet-cw0)*                   &
          dsin(kw*(xwu-cw0*Time))*dcosh(kw*(PGrid%y(1,j)-HChannel))/dsinh(kw*Ha)
        end if
        if(PGrid%y(1,j)-Hw<ywv) then
          Vari%v(ibeg-ight,j)=Amp0*kw*(UwInlet-cw0)*dcos(kw*(xwv-cw0*Time))*   &
                  dsinh(kw*(PGrid%y(1,j)+0.5d0*PGrid%dy(1,j)))/dsinh(kw*Hw)
        else
          Vari%v(ibeg-ight,j)=-Amp0*kw*(UgInlet-cw0)*dcos(kw*(xwv-cw0*Time))*  &
              dsinh(kw*(PGrid%y(1,j)+0.5d0*PGrid%dy(1,j)-HChannel))/dsinh(kw*Ha)
        end if
      ! Right outlet water
        Vari%p(Isize+ibeg+ight-1,j)=0.d0!Vari%p(Isize+ibeg-1,j)
        Vari%t(Isize+ibeg+ight-1,j)=Vari%t(ibeg+Isize-1,j)
        Vari%u(Isize+ibeg+ight-1,j)=Vari%u(Isize+ibeg-1,j) !Vari%Uint/Vari%Uref
        Vari%v(Isize+ibeg+ight-1,j)=Vari%v(Isize+ibeg-1,j) !Vari%Vint/Vari%Uref
      end do
      do i = ibeg,Isize+ibeg-1
      ! Bottom wall-slip wall
        Vari%p(i,jbeg-jght)=Vari%p(i,jbeg)
        Vari%t(i,jbeg-jght)=Vari%t(i,jbeg)
        Vari%u(i,jbeg-jght)=Vari%u(i,jbeg)
        Vari%v(i,jbeg-jght)=0.d0 !Vari%Vint/Vari%Uref
      ! Open air
        Vari%p(i,Jsize+jbeg+jght-1)=0.d0
        Vari%t(i,jbeg+Jsize+jght-1)=Vari%t(i,jbeg+Jsize-1)
        Vari%u(i,Jsize+jbeg+jght-1)=Vari%u(i,Jsize+jbeg-1)
     !   Vari%v(i,Jsize+jbeg-1) = Vari%v(i,Jsize+jbeg-2)
        Vari%v(i,Jsize+jght) = Vari%v(i,Jsize)
      end do
    end subroutine Boundary_Condition_Var
End Module StateVariables



