Module Particles
  USE PrecisionVar
  USE MPI
  USE Mesh, ONLY : TPoint, Grid, Cell, ibeg, jbeg, Isize, Jsize, ight, jght
  USE StateVariables, ONLY : TVariables, TWave
  USE Constants, ONLY : pi, epsi, nua, nuw, roa, row
  USE Clsvof,ONLY: vofeps,SolidObject
  use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
  IMPLICIT NONE
  PRIVATE
  INTEGER(kind=it4b),PARAMETER:: itp=10
  REAL(KIND=dp),PARAMETER:: tol=1.d-20,Ca=0.5d0
  TYPE,PRIVATE:: ParVel
    REAL(KIND=dp):: u,v
  End TYPE ParVel
  TYPE :: TParticle
    INTEGER(kind=it4b):: np
    INTEGER(it4b) :: NParInlet, IParInlet
    TYPE(TPoint),DIMENSION(:),allocatable:: Posp
    TYPE(ParVel),DIMENSION(:),allocatable:: uvp
    REAL(dp),DIMENSION(:),allocatable:: mp,dp,tp
    REAL(dp) :: UParInlet, HParInlet, DParInlet, zp, rop, gx, gy
  ! For secondary break up
    REAL(KIND=dp),DIMENSION(:),allocatable:: t
    REAL(KIND=dp),DIMENSION(:),allocatable:: VrelG,y,dydt
  CONTAINS
    PROCEDURE, PASS(this), PUBLIC :: InitializeParticles
    PROCEDURE, PASS(this), PUBLIC :: ParticleInletCondition
    PROCEDURE, PASS(this), PUBLIC :: TrackingParticles
  End TYPE TParticle
  PUBLIC:: ParticlePosition,Drag
  Interface ParticlePosition
    Module procedure ParticlePosition
  End interface
  Interface Drag
    Module procedure Drag
  End interface
  Interface TParticle
     Module procedure construct
  End interface TParticle

  PUBLIC :: TParticle
  Contains

  TYPE(TParticle) function construct(np, NParInlet, IParInlet) RESULT ( this )
    !
    INTEGER(it4b), intent(in) :: np, NParInlet, IParInlet
    INTEGER (it4b) :: np_
    this % np = np
    this % NParInlet = NParInlet
    this % IParInlet = IParInlet
    np_ = MAX(1, np+NParInlet)
    allocate(this%dp(np_))
    allocate(this%Posp(np_))
    allocate(this%uvp(np_))
    allocate(this%mp(np_))
    allocate(this%tp(np_))
    allocate(this%t(np_))
    allocate(this%VrelG(np_))
    allocate(this%y(np_))
    allocate(this%dydt(np_))
  end function construct
  Subroutine InitializeParticles(this, Pgrid,Var, UParInlet, HParInlet, DParInlet, zp, rop, gx, gy)
    IMPLICIT NONE
    CLASS(TParticle),INTENT(INOUT):: this
    TYPE(Grid),INTENT(IN):: PGrid
    TYPE(TVariables),INTENT(IN):: Var
    REAL(dp), INTENT(in) :: UParInlet, HParInlet, DParInlet, zp, rop, gx, gy
    INTEGER:: i
    REAL(KIND=dp):: DragFC
    REAL(KIND=dp),DIMENSION(:),allocatable:: ranum
    allocate(ranum(this%np))
    call Random_Number(ranum)
    DragFC=18.d0
    this % UParInlet = UParInlet
    this % HParInlet = HParInlet
    this % DParInlet = DParInlet
    this % zp        = zp
    this % rop       = rop
    this % gx        = gx
    this % gy        = gy
    do i=1,this%np
      this%dp(i)=1.d-6+0.008d0*ranum(i)
    end do
    call Random_Number(ranum)
    do i=1,this%np
      this%Posp(i)%x=(PGrid%x(10,1)+(PGrid%x(Isize,1)-PGrid%x(350,1))*ranum(i))*Pgrid%Lref
    end do
    call Random_Number(ranum)
    do i=1,this%np
      this%Posp(i)%y=(PGrid%y(1,2)+(PGrid%y(1,Jsize)-PGrid%y(1,80))*ranum(i))*Pgrid%Lref
    end do
 !  No contribution of particles
 !******************************
    do i=1,this%np
      this%mp(i)=1.d0/6.d0*pi*(this%dp(i))**3.d0*row
      this%tp(i)=row*this%dp(i)**2.d0/DragFC/(nua*roa)
      this%uvp(i)%u=20.d0
      this%uvp(i)%v=0.d0
    end do
    deallocate(ranum)
  end subroutine InitializeParticles
  subroutine TrackingParticles(this, PGrid,PCell,Var,dt,BoomCase)
    IMPLICIT NONE
    CLASS(TParticle),INTENT(INOUT):: this 
    TYPE(Grid),INTENT(IN):: PGrid
    TYPE(Cell),INTENT(IN):: PCell
    TYPE(TVariables),INTENT(IN):: Var
    REAL(KIND=dp),INTENT(IN):: dt
    TYPE(SolidObject),INTENT(IN),optional :: BoomCase
    INTEGER(kind=it4b):: i,ii,jj,itep
    REAL(KIND=dp),DIMENSION(:),allocatable:: Upo,Upn,Vpo,Vpn,axp,ayp
    TYPE(TPoint),DIMENSION(:),allocatable:: xyp
    REAL(KIND=dp):: dudx,dvdy,ug,vg,ug0,vg0,VRel,Reyp,Cd
    REAL(KIND=dp):: tp,FXT,EXPT,Spx,Spy,dtp,gama,nupp,ropp,beta,sig
    allocate(Upo(this%np))
    allocate(Vpo(this%np))
    allocate(Upn(this%np))
    allocate(Vpn(this%np))
    allocate(xyp(this%np))
    allocate(axp(this%np))
    allocate(ayp(this%np))
    dtp = dt*(PGrid%Lref/Var%Uref)/dble(itp)
    do i=1,this%np
      xyp(i)%x=this%Posp(i)%x
      xyp(i)%y=this%PosP(i)%y
      Upo(i)=this%uvp(i)%u
      VPo(i)=this%uvp(i)%v
    end do
    do itep=1,itp
      do i=1,this%np
        call ParticlePosition(this%Posp(i),PGrid,ii,jj)
        if(ii/=-1.and.jj/=-1) then
          if(this%PosP(i)%x/Pgrid%Lref<PGrid%x(Isize,1).and.                       &
             this%PosP(i)%y/Pgrid%Lref<PGrid%y(1,Jsize).and.                       &
             PCell%vofS(ii,jj)<1.d0-epsi) then
            ropp=row*PCell%vof(ii,jj)/(1.d0-PCell%vofS(ii,jj)+tol)+            &
                 roa*(1.d0-PCell%vof(ii,jj)-PCell%vofS(ii,jj))/                &
                     (1.d0-PCell%vofS(ii,jj)+tol)
            gama=this%rop/ropp
            nupp=nuw*PCell%vof(ii,jj)/(1.d0-PCell%vofS(ii,jj)+tol)+            &
                 nua*(1.d0-PCell%vof(ii,jj)-PCell%vofS(ii,jj))/                &
                     (1.d0-PCell%vofS(ii,jj)+tol)
            dudx=(Var%u(ii,jj)-Var%u(ii-1,jj))/PGrid%dx(ii,jj)
            dvdy=(Var%v(ii,jj)-Var%v(ii,jj-1))/PGrid%dy(ii,jj)
            ug0=0.5d0*(Var%u(ii,jj)+Var%u(ii-1,jj))
            vg0=0.5d0*(Var%v(ii,jj)+Var%v(ii,jj-1))
            ug=(ug0+dudx*(this%Posp(i)%x/Pgrid%Lref-PGrid%x(ii,jj)))*Var%Uref
            vg=(vg0+dvdy*(this%Posp(i)%y/Pgrid%Lref-PGrid%y(ii,jj)))*Var%Uref
            VRel=dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
            if(VRel<1.d-7) VRel=1.d-7
            Reyp=this%dp(i)*VRel/nupp
            Cd=Drag(Reyp)
            axp(i)=(1.d0-gama)/(gama+Ca)*this%gx-3.d0*Cd/8.d0/(this%dp(i)/2.d0)/   &
                   (gama+Ca)*(upo(i)-ug)*dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
            ayp(i)=(1.d0-gama)/(gama+Ca)*this%gy-3.d0*Cd/8.d0/(this%dp(i)/2.d0)/   &
                (gama+Ca)*(vpo(i)-vg)*dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
         !   Upn(i)=Upo(i)+axp(i)*dtp
         !   Vpn(i)=Vpo(i)+ayp(i)*dtp
            sig=(1.d0-gama)/(gama+Ca)*this%gx
            beta=3.d0*Cd/8.d0/(this%dp(i)/2.d0)/(gama+Ca)*                    &
                                      dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
            Upn(i)=Upo(i)*exp(-beta*dtp)+ug*(1.d0-exp(-beta*dtp))+              &
                   sig/beta*(1.d0-exp(-beta*dtp))
            xyp(i)%x=xyp(i)%x+0.5d0*(Upn(i)+Upo(i))*dtp
            sig=(1.d0-gama)/(gama+Ca)*this%gy
            Vpn(i)=Vpo(i)*exp(-beta*dtp)+vg*(1.d0-exp(-beta*dtp))+              &
                   sig/beta*(1.d0-exp(-beta*dtp))
            xyp(i)%y=xyp(i)%y+0.5d0*(Vpn(i)+Vpo(i))*dtp
            if(isnan(Upo(i)).or.dabs(xyp(i)%x)>1.d10.or.dabs(axp(i))>1.d10.or.  &
                                            isnan(upn(i)).or.isnan(vpn(i))) then
              print*,i
              print*,Cd,beta,gama
              
              print*,Reyp, nupp, Vrel 
              print*, 
              print*,axp(i),upn(i),vpn(i)
              print*,
              print*,axp(i),dtp
              print*,(1.d0-gama)/(gama+Ca)*this%gx,3.d0*Cd/8.d0/(this%dp(i)/2.d0)/   &
                   (gama+Ca)*(upo(i)-ug)*dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
              print*,3.d0*Cd/8.d0,(this%dp(i)/2.d0)
              print*, '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
              print*,(gama+Ca)*(upo(i)-ug)*dsqrt((upo(i)-ug)**2.d0+(vpo(i)-vg)**2.d0)
              print*,Reyp
              print*,'What the fuck'
              print*,this%uvp(i)%u,this%uvp(i)%v
              print*,nupp
              print*,'-----------------------------------------------------------'
              print*,PCell%vof(ii,jj)/(1.d0-PCell%vofS(ii,jj)+tol)
              print*,(1.d0-PCell%vof(ii,jj)-PCell%vofS(ii,jj))/                &
                     (1.d0-PCell%vofS(ii,jj)+tol)
              print*,Upo(i),exp(-beta*dtp),ug*(1.d0-exp(-beta*dtp))
              print*,sig,beta,(1.d0-exp(-beta*dtp))
              print*,upn(i),upo(i)
              print*, '111111111111111111111111111111111111111111111+++++++++++++'
              print*,dudx,dvdy
              print*,beta,sig,exp(-beta*dtp)
              print*,(1.d0-gama),(gama+Ca)
              print*,gama,this%gx
              print*,Cd,(this%dp(i)/2.d0)
              print*,upo(i),ug
              print*,ug0,Var%u(ii,jj),Var%u(ii-1,jj)
              print*,this%Posp(i)%x,this%Posp(i)%y
              print*,PCell%vof(ii,jj)
              print*, 'Particle 127'
              read(stdin,*)
            end if
            Upo(i)=Upn(i)
            Vpo(i)=Vpn(i)
          elseif(PCell%vofS(ii,jj)>1.d0-epsi) then
            if(present(BoomCase)) then 
              Upn(i)=BoomCase%us
              Vpn(i)=BoomCase%vs
            end if
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
    do i=1,this%np
      this%Posp(i)%x=xyp(i)%x
      this%PosP(i)%y=xyp(i)%y
      this%uvp(i)%u=Upo(i)
      this%uvp(i)%v=VPo(i)
      if(isnan(this%uvp(i)%u)) then
        print*, 'Particle 149'
        read(stdin,*)
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
    TYPE(TPoint),INTENT(IN)::Posi
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

  Subroutine ParticleInletCondition(this, PGrid,PCell,wave)
    IMPLICIT NONE
    CLASS(TParticle),INTENT(INOUT):: this
    TYPE(Grid),INTENT(IN):: PGrid
    TYPE(Cell),INTENT(IN):: PCell
    TYPE(TWave), INTENT(in) :: wave
    INTEGER(KIND=it4b):: i,j
    REAL(KIND=dp),DIMENSION(:),allocatable:: ranum
    allocate(ranum(this%NParInlet))
    call Random_Number(ranum)
    do i=this%np+1,this%np+this%NParInlet
      this%dp(i) = this % DParInlet*(ranum(i-this%np)+1.d-3)
    end do
    call Random_Number(ranum)
    do i=this%np+1,this%np+this%NParInlet
      this%Posp(i)%x=PGrid%x(1,1)+PGrid%dx(1,1)/3.d0
      this%Posp(i)%y=wave%Depthw-wave%Amp0*2.d0-this%HParInlet*ranum(i-this%np)
    end do
    call Random_Number(ranum)
    do i=this%np+1,this%np+this%NParInlet
      this%uvp(i)%u = this % UParInlet*(ranum(i-this%np)+0.5d0)
      this%uvp(i)%v = 0.d0
    end do
    this%np=this%np+this%NParInlet !FIXME: This does not make sense anymore, revise this and initialize
  End subroutine ParticleInletCondition
End module Particles
