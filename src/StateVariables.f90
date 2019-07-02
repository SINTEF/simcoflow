Module StateVariables
    USE PrecisionVar
    USE Mesh
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
    LOGICAL,PUBLIC:: RunAgain,ICorProb
!   For particles from 3D to 2D
    REAL(KIND=dp),PUBLIC::zp,UParInlet,HParInlet,DParInlet,rop,gx,gy
!   For number of particles into system
    INTEGER(kind=it4b),PUBLIC:: NParInlet,IParInlet
!   For Water Entry problem
    REAL(KIND=dp),PUBLIC::y0,y1
!   For Wave Braking on vertical wall
    REAL(KIND=dp),PUBLIC::t0,cw0,Amp0,Depthw,Lamdaw,twp,HChannel,              &
                          LChannel,kw,omew,HDomain
    REAL(KIND=dp),PUBLIC,PARAMETER::nuw=1.0034d-6,nua=1.506d-5,                &
                                    roa=1.225d0,row=998.2d0 ! At T = 20oC
    REAL(KIND=dp),PUBLIC,PARAMETER::epsi=1.d-3,epsiF=1.d-2,BetaVis=0.5d0
    CHARACTER*70,PUBLIC           ::dir
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
          if(PGrid%y(i,j)-Depthw<ywu) then
            Vari%u(i,j)=UwInlet-Amp0*kw*(UwInlet-cw0)*dsin(kw*xwu)*            &
                        dcosh(kw*PGrid%y(i,j))/dsinh(kw*Depthw)
          else
            Vari%u(i,j)=UgInlet+Amp0*kw*(UgInlet-cw0)*dsin(kw*xwu)*            &
                               dcosh(kw*(PGrid%y(i,j)-HChannel))/dsinh(kw*Ha)
          end if
          if(PGrid%y(i,j)-Depthw<ywv) then
            Vari%v(i,j)=Amp0*kw*(UwInlet-cw0)*dcos(kw*xwv)*                    &
               dsinh(kw*(PGrid%y(i,j)+0.5d0*PGrid%dy(i,j)))/dsinh(kw*Depthw)
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
   !   call Boundary_Condition_Var(PGrid,PCell,Vari,0.d0)
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
    subroutine Boundary_Condition_Var(PGrid,PCell,Vari,Time)
      TYPE(Grid),INTENT(IN):: PGrid
      TYPE(Cell),INTENT(IN):: PCell
      TYPE(Variables),INTENT(INOUT):: Vari
      REAL(KIND=dp),INTENT(IN):: Time
      INTEGER(kind=it4b):: i,j,temp
      REAL(KIND=dp),PARAMETER:: Twall = 300.d0
      REAL(KIND=dp):: ywu,ywv,xwu,xwv,ywuout,ywvout,xwuout,xwvout,            &
                      Hwin,Hwout,Hain,Haout
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
      ywu=Amp0*dsin(kw*(xwu-cw0*Time))
      ywv=Amp0*dsin(kw*(xwv-cw0*Time))

      xwuout=PGrid%x(Isize,1)+PGrid%dx(Isize,1)/2.d0
      xwvout=PGrid%x(Isize,1)+PGrid%dx(Isize,1)
      ywuout=Amp0*dsin(kw*(xwuout-cw0*Time))
      ywvout=Amp0*dsin(kw*(xwvout-cw0*Time))

      if(temp==2) Hwout=Depthw
      Hain=HChannel-Hwin
      Haout=HChannel-Hwout
      do j = jbeg,Jsize+jbeg-1
      ! Left inlet water
        Vari%p(ibeg-ight,j)=Vari%p(ibeg,j) !Vari%Pint/(Vari%Pref)
        Vari%t(ibeg-ight,j)=Vari%t(ibeg,j)
        if(PGrid%y(1,j)-Hwin<ywu) then
          Vari%u(ibeg-ight,j)=UwInlet-Amp0*kw*(UwInlet-cw0)*                   &
                    dsin(kw*(xwu-cw0*Time))*dcosh(kw*PGrid%y(1,j))/dsinh(kw*Hwin)
        else
          Vari%u(ibeg-ight,j)=UgInlet+Amp0*kw*(UgInlet-cw0)*                   &
          dsin(kw*(xwu-cw0*Time))*dcosh(kw*(PGrid%y(1,j)-HChannel))/dsinh(kw*Hain)
        end if

        if(PGrid%y(1,j)-Hwin<ywv) then
          Vari%v(ibeg-ight,j)=Amp0*kw*(UwInlet-cw0)*dcos(kw*(xwv-cw0*Time))*   &
                  dsinh(kw*(PGrid%y(1,j)+0.5d0*PGrid%dy(1,j)))/dsinh(kw*Hwin)
        else
          Vari%v(ibeg-ight,j)=-Amp0*kw*(UgInlet-cw0)*dcos(kw*(xwv-cw0*Time))*  &
              dsinh(kw*(PGrid%y(1,j)+0.5d0*PGrid%dy(1,j)-HChannel))/dsinh(kw*Hain)
        end if

      ! Right outlet water
      !  if(PGrid%y(Isize,j)-Hwout<ywuout) then
        if(PCell%vof(isize,j)>0.5d0) then
          Vari%u(Isize,j)=UwInlet!-Amp0*kw*(UwInlet-cw0)*                       &
            ! dsin(kw*(xwuout-cw0*Time))*dcosh(kw*PGrid%y(Isize,j))/dsinh(kw*Hwout)
        else
          Vari%u(Isize,j)=UgInlet  !+Amp0*kw*(UgInlet-cw0)*                       &
            ! dsin(kw*(xwuout-cw0*Time))*dcosh(kw*(PGrid%y(Isize,j)-HChannel))/ &
            !                                                     dsinh(kw*Haout)
        end if

        if(PGrid%y(Isize,j)-Hwout<ywvout) then
          Vari%v(Isize+1,j)=Amp0*kw*(UwInlet-cw0)*dcos(kw*(xwvout-cw0*Time))*  &
              dsinh(kw*(PGrid%y(Isize,j)+0.5d0*PGrid%dy(Isize,j)))/dsinh(kw*Hwout)
        else
          Vari%v(Isize+1,j)=-Amp0*kw*(UgInlet-cw0)*dcos(kw*(xwvout-cw0*Time))* &
              dsinh(kw*(PGrid%y(Isize,j)+0.5d0*PGrid%dy(Isize,j)-HChannel))/   &
                                                                  dsinh(kw*Haout)
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
        Vari%u(i,jbeg-jght)=-Amp0*kw*(UwInlet-cw0)*dsin(kw*(xwu-cw0*Time))*    &
                        dcosh(kw*(PGrid%y(i,1)-PGrid%dy(i,1)))/dsinh(kw*Depthw)
        xwv=PGrid%x(i,1)
        Vari%v(i,jbeg-jght)=Amp0*kw*(UwInlet-cw0)*dcos(kw*(xwv-cw0*Time))*    &
                  dsinh(kw*(PGrid%y(i,1)-0.5d0*PGrid%dy(i,1)))/dsinh(kw*Depthw) !0.d0 !Vari%Vint/Vari%Uref
      ! Open air
        Vari%p(i,Jsize+jbeg+jght-1)=0.d0
        Vari%t(i,jbeg+Jsize+jght-1)=Vari%t(i,jbeg+Jsize-1)
        Vari%u(i,Jsize+jbeg+jght-1)=Vari%u(i,Jsize+jbeg-1)
     !   Vari%v(i,Jsize+jbeg-1) = Vari%v(i,Jsize+jbeg-2)
        Vari%v(i,Jsize+jght) = Vari%v(i,Jsize)
      end do
    end subroutine Boundary_Condition_Var
End Module StateVariables



