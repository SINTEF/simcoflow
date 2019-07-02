Module Particles
  USE PrecisionVar
  USE MPI
  USE Mesh
  USE StateVariables
  USE Clsvof,ONLY: vofeps,SolidObject
  PRIVATE
  INTEGER(kind=it4b),PARAMETER:: itp=10
  REAL(KIND=dp),PARAMETER:: tol=1.d-20,Ca=0.5d0
  TYPE,PRIVATE:: ParVel
    REAL(KIND=dp):: u,v
  End TYPE ParVel
  TYPE,PUBLIC:: Particle
    INTEGER(kind=it4b):: np
    TYPE(Point),DIMENSION(:),allocatable:: Posp
    TYPE(ParVel),DIMENSION(:),allocatable:: uvp
    REAL(KIND=dp),DIMENSION(:),allocatable:: mp,dp,tp
  ! For secondary break up
    REAL(KIND=dp),DIMENSION(:),allocatable:: t
    REAL(KIND=dp),DIMENSION(:),allocatable:: VrelG,y,dydt
  End TYPE Particle
  PUBLIC:: TrackingParticles,InitializeParticles,ParticlePosition,Drag,ParticleInletCondition
  Interface TrackingParticles
    Module Procedure TrackingParticles
  End interface
  Interface InitializeParticels
    Module procedure InitializeParticles
  End interface
  Interface ParticlePosition
    Module procedure ParticlePosition
  End interface
  Interface ParticleInletCondition
    Module procedure ParticleInletCondition
  End interface
  Interface Drag
    Module procedure Drag
  End interface
  Contains
  Subroutine InitializeParticles(Pgrid,Var,TraPar)
    IMPLICIT NONE
    TYPE(Grid),INTENT(IN):: PGrid
    TYPE(Variables),INTENT(IN):: Var
    TYPE(Particle),INTENT(INOUT):: TraPar
    INTEGER:: i
    REAL(KIND=dp):: DragFC
    REAL(KIND=dp),DIMENSION(:),allocatable:: ranum
    allocate(ranum(TraPar%np))
    call Random_Number(ranum)
    DragFC=18.d0
    do i=1,TraPar%np
      TraPar%dp(i)=1.d-6+0.008d0*ranum(i)
    end do
    call Random_Number(ranum)
    do i=1,TraPar%np
      TraPar%Posp(i)%x=(PGrid%x(10,1)+(PGrid%x(Isize,1)-PGrid%x(350,1))*ranum(i))*Lref
    end do
    call Random_Number(ranum)
    do i=1,TraPar%np
      TraPar%Posp(i)%y=(PGrid%y(1,2)+(PGrid%y(1,Jsize)-PGrid%y(1,80))*ranum(i))*Lref
    end do
 !  No contribution of particles
 !******************************
    do i=1,TraPar%np
      TraPar%mp(i)=1.d0/6.d0*pi*(TraPar%dp(i))**3.d0*row
      TraPar%tp(i)=row*TraPar%dp(i)**2.d0/DragFC/(nua*roa)
      TraPar%uvp(i)%u=20.d0
      TraPar%uvp(i)%v=0.d0
    end do
    deallocate(ranum)
  end subroutine InitializeParticles
  subroutine TrackingParticles(PGrid,PCell,BoomCase,Var,dt,TraPar)
    IMPLICIT NONE
    TYPE(Grid),INTENT(IN):: PGrid
    TYPE(Cell),INTENT(IN):: PCell
    TYPE(SolidObject),INTENT(INOUT):: BoomCase
    TYPE(Variables),INTENT(IN):: Var
    REAL(KIND=dp),INTENT(IN):: dt
    TYPE(Particle),INTENT(INOUT):: TraPar
    INTEGER(kind=it4b):: i,ii,jj,itep
    REAL(KIND=dp),DIMENSION(:),allocatable:: Upo,Upn,Vpo,Vpn,axp,ayp
    TYPE(Point),DIMENSION(:),allocatable:: xyp
    REAL(KIND=dp):: dudx,dvdy,ug,vg,ug0,vg0,VRel,Reyp,Cd
    REAL(KIND=dp):: tp,FXT,EXPT,Spx,Spy,dtp,gama,nupp,ropp,beta,sig
    allocate(Upo(TraPar%np))
    allocate(Vpo(TraPar%np))
    allocate(Upn(TraPar%np))
    allocate(Vpn(TraPar%np))
    allocate(xyp(TraPar%np))
    allocate(axp(TraPar%np))
    allocate(ayp(TraPar%np))
    dtp = dt*(PGrid%Lref/Var%Uref)/dble(itp)
    do i=1,TraPar%np
      xyp(i)%x=TraPar%Posp(i)%x
      xyp(i)%y=TraPar%PosP(i)%y
      Upo(i)=TraPar%uvp(i)%u
      VPo(i)=TraPar%uvp(i)%v
    end do
    do itep=1,itp
      do i=1,TraPar%np
        call ParticlePosition(TraPar%Posp(i),PGrid,ii,jj)
        if(ii/=-1.and.jj/=-1) then
          if(TraPar%PosP(i)%x/Lref<PGrid%x(Isize,1).and.                       &
             TraPar%PosP(i)%y/Lref<PGrid%y(1,Jsize).and.                       &
             PCell%vofS(ii,jj)<1.d0-epsi) then
            ropp=row*PCell%vof(ii,jj)/(1.d0-PCell%vofS(ii,jj)+tol)+            &
                 roa*(1.d0-PCell%vof(ii,jj)-PCell%vofS(ii,jj))/                &
                     (1.d0-PCell%vofS(ii,jj)+tol)
            gama=rop/ropp
            nupp=nuw*PCell%vof(ii,jj)/(1.d0-PCell%vofS(ii,jj)+tol)+            &
                 nua*(1.d0-PCell%vof(ii,jj)-PCell%vofS(ii,jj))/                &
                     (1.d0-PCell%vofS(ii,jj)+tol)
            dudx=(Var%u(ii,jj)-Var%u(ii-1,jj))/PGrid%dx(ii,jj)
            dvdy=(Var%v(ii,jj)-Var%v(ii,jj-1))/PGrid%dy(ii,jj)
            ug0=0.5d0*(Var%u(ii,jj)+Var%u(ii-1,jj))
            vg0=0.5d0*(Var%v(ii,jj)+Var%v(ii,jj-1))
            ug=(ug0+dudx*(TraPar%Posp(i)%x/Lref-PGrid%x(ii,jj)))*Var%Uref
            vg=(vg0+dvdy*(TraPar%Posp(i)%y/Lref-PGrid%y(ii,jj)))*Var%Uref
            VRel=dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
            if(VRel<1.d-7) VRel=1.d-7
            Reyp=TraPar%dp(i)*VRel/nupp
            Cd=Drag(Reyp)
            axp(i)=(1.d0-gama)/(gama+Ca)*gx-3.d0*Cd/8.d0/(TraPar%dp(i)/2.d0)/   &
                   (gama+Ca)*(upo(i)-ug)*dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
            ayp(i)=(1.d0-gama)/(gama+Ca)*gy-3.d0*Cd/8.d0/(TraPar%dp(i)/2.d0)/   &
                (gama+Ca)*(vpo(i)-vg)*dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
         !   Upn(i)=Upo(i)+axp(i)*dtp
         !   Vpn(i)=Vpo(i)+ayp(i)*dtp
            sig=(1.d0-gama)/(gama+Ca)*gx
            beta=3.d0*Cd/8.d0/(TraPar%dp(i)/2.d0)/(gama+Ca)*                    &
                                      dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
            Upn(i)=Upo(i)*exp(-beta*dtp)+ug*(1.d0-exp(-beta*dtp))+              &
                   sig/beta*(1.d0-exp(-beta*dtp))
            xyp(i)%x=xyp(i)%x+0.5d0*(Upn(i)+Upo(i))*dtp
            sig=(1.d0-gama)/(gama+Ca)*gy
            Vpn(i)=Vpo(i)*exp(-beta*dtp)+vg*(1.d0-exp(-beta*dtp))+              &
                   sig/beta*(1.d0-exp(-beta*dtp))
            xyp(i)%y=xyp(i)%y+0.5d0*(Vpn(i)+Vpo(i))*dtp
            if(isnan(Upo(i)).or.dabs(xyp(i)%x)>1.d10.or.dabs(axp(i))>1.d10.or.  &
                                            isnan(upn(i)).or.isnan(vpn(i))) then
              print*,i
              print*,axp(i),dtp
              print*,
              print*,(1.d0-gama)/(gama+Ca)*gx,3.d0*Cd/8.d0/(TraPar%dp(i)/2.d0)/   &
                   (gama+Ca)*(upo(i)-ug)*dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
              print*,3.d0*Cd/8.d0,(TraPar%dp(i)/2.d0)
              print*,(gama+Ca)*(upo(i)-ug)*dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
              print*,
              print*,Reyp
              print*,
              print*,'What the fuck'
              print*,TraPar%uvp(i)%u,TraPar%uvp(i)%v
              print*,nupp
              print*,PCell%vof(ii,jj)/(1.d0-PCell%vofS(ii,jj)+tol)
              print*,(1.d0-PCell%vof(ii,jj)-PCell%vofS(ii,jj))/                &
                     (1.d0-PCell%vofS(ii,jj)+tol)
              print*,
              print*,Upo(i),exp(-beta*dtp),ug*(1.d0-exp(-beta*dtp))
              print*,sig,beta,(1.d0-exp(-beta*dtp))
              print*,upn(i),upo(i)
              print*,
              print*,dudx,dvdy
              print*,
              print*,beta,sig,exp(-beta*dtp)
              print*,
              print*,(1.d0-gama),(gama+Ca)
              print*,gama,gx
              print*,Cd,(TraPar%dp(i)/2.d0)
              print*,upo(i),ug
              print*,ug0,Var%u(ii,jj),Var%u(ii-1,jj)
              print*,
              print*,TraPar%Posp(i)%x,TraPar%Posp(i)%y
              print*,PCell%vof(ii,jj)
              pause 'Particle 127'
            end if
            Upo(i)=Upn(i)
            Vpo(i)=Vpn(i)
          elseif(PCell%vofS(ii,jj)>1.d0-epsi) then
            Upn(i)=BoomCase%us
            Vpn(i)=BoomCase%vs
            xyp(i)%x=xyp(i)%x+0.5d0*(Upn(i)+Upo(i))*dtp
            xyp(i)%y=xyp(i)%y+0.5d0*(Vpn(i)+Vpo(i))*dtp
            Upo(i)=Upn(i)
            Vpo(i)=Vpn(i)
          else
            Upn(i)=0.d0
            Vpn(i)=0.d0
            xyp(i)%x=1.d3
            xyp(i)%y=1.d3
          end if
        end if
      end do
    end do
    do i=1,TraPar%np
      TraPar%Posp(i)%x=xyp(i)%x
      TraPar%PosP(i)%y=xyp(i)%y
      TraPar%uvp(i)%u=Upo(i)
      TraPar%uvp(i)%v=VPo(i)
      if(isnan(TRaPar%uvp(i)%u)) then
        pause 'Particle 149'
      end if
    end do
    deallocate(Upo,Upn)
    deallocate(Vpo,Vpn)
    deallocate(xyp,axp,ayp)
  end Subroutine TrackingParticles
  Function Drag(Reyp) result(Cd)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!                                                              !
!    'DRAG' CALCULATES DRAG COEFFICIENT                        !
!                                                              !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
    REAL(KIND=dp):: Cd,Reyp
    if(Reyp<1.d0) then
      Cd = 24.d0/(Reyp+tol)
      return
    end if
    if(Reyp<=1000.d0) then
      Cd = 24.d0/(Reyp+tol)*(1.d0+0.15d0*Reyp**0.687d0)
      return
    end if
    Cd = 0.44d0
    return
  end function Drag
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
! Calculate the position of particle                            !
! The output is the cell number                                 !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
  Subroutine ParticlePosition(Posi,TGrid,ii,jj)
    TYPE(Point),INTENT(IN)::Posi
    TYPE(Grid),INTENT(IN)::TGrid
    INTEGER(kind=it4b),INTENT(OUT):: ii,jj
    INTEGER(kind=it4b):: i,j
    ii=-1
    do i=1,Isize
      if(Posi%x/TGrid%Lref>=TGrid%x(i,1)-TGrid%dx(i,1)/2.d0.and.               &
         Posi%x/TGrid%Lref<=TGrid%x(i,1)+TGrid%dx(i,1)/2.d0) then
        ii=i
        exit
      end if
    end do
    jj=-1
    do j=1,Jsize
      if(Posi%y/TGrid%Lref>=TGrid%y(1,j)-TGrid%dy(1,j)/2.d0.and.               &
         Posi%y/TGrid%Lref<=TGrid%y(1,j)+TGrid%dy(1,j)/2.d0) then
        jj=j
        exit
      end if
    end do
  End subroutine ParticlePosition

  Subroutine ParticleInletCondition(PGrid,PCell,TraPar)
    IMPLICIT NONE
    TYPE(Grid),INTENT(IN):: PGrid
    TYPE(Cell),INTENT(IN):: PCell
    TYPE(Particle),INTENT(INOUT):: TraPar
    INTEGER(KIND=it4b):: i,j
    REAL(KIND=dp),DIMENSION(:),allocatable:: ranum
    allocate(ranum(NParInlet))
    call Random_Number(ranum)
    do i=TraPar%np+1,TraPar%np+NParInlet
      TraPar%dp(i)=DParInlet*(ranum(i-TraPar%np)+1.d-3)
    end do
    call Random_Number(ranum)
    do i=TraPar%np+1,TraPar%np+NParInlet
      TraPar%Posp(i)%x=PGrid%x(1,1)+PGrid%dx(1,1)/3.d0
      TraPar%Posp(i)%y=Depthw-Amp0*2.d0-HParInlet*ranum(i-TraPar%np)
    end do
    call Random_Number(ranum)
    do i=TraPar%np+1,TraPar%np+NParInlet
      TraPar%uvp(i)%u=UParInlet*(ranum(i-TraPar%np)+0.5d0)
      TraPar%uvp(i)%v=0.d0
    end do
    TraPar%np=TraPar%np+NParInlet
  End subroutine ParticleInletCondition
End module Particles
