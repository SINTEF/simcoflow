Module PredictorUV
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE Clsvof
    USE StateVariables
    USE Matrix
    USE Printresult
    USE MPI
    USE Particles
    IMPLICIT NONE
    PRIVATE
    real(dp),DIMENSION(:,:),pointer:: u,Uvolf
    real(dp),DIMENSION(:,:),pointer:: v,Vvolf
    Real(dp),DIMENSION(:,:),pointer:: p
    Real(dp),PARAMETER,PRIVATE:: tol=1.d-14
    TYPE,PUBLIC:: Predictor
        REAL(KIND=dp),DIMENSION(:,:),allocatable:: u,v
    End TYPE
    TYPE,PUBLIC:: PoissonCoefficient
        REAL(KIND=dp),DIMENSION(:,:),allocatable:: Dp
    End TYPE
    PUBLIC:: Predictor_UV
    Interface Predictor_UV
        Module Procedure Predictor_UV
    End interface
    Contains
    Subroutine Predictor_UV(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,PCell,      &
                            UCell,VCell,TVar,PU,PV,Pred,Flux_n1,TraPar,        &
                            VolParU,VolParV,SParU,SParV,BoomCase,dt,itt)
      IMPLICIT NONE
      INTEGER(kind=it8b),INTENT(IN):: itt
      REAL(KIND=dp),INTENT(IN):: dt
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: PCell,PCellO
      TYPE(Cell),INTENT(IN),target:: UCell,VCell,UCellO,VCellO
      TYPE(Variables),INTENT(IN),target:: TVar
      TYPE(Predictor),INTENT(INOUT):: Pred
      TYPE(PoissonCoefficient),INTENT(INOUT):: PU,PV
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT):: Flux_n1
      REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(IN):: VolParU,VolParV,SParU,SParV
      TYPE(Particle),INTENT(INOUT):: TraPar
      TYPE(SolidObject),INTENT(IN):: BoomCase
      INTEGER(kind=it4b):: i,j,ii,jj
      INTEGER*8:: A,parcsr_A,b,par_b,x,par_x,solver,precond
      INTEGER(kind=it4b),DIMENSION(:,:),allocatable:: MoParExIJU,MoParExIJV
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: ConFluxEW,ConFluxNS,        &
                    DifFluxEW,DifFluxNS,ExEDFluxEW,ExEDFluxNS,EDFluxEW,EDFluxNS
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: CFluxEW,CFluxNS,            &
                    MaFluxEW,MaFluxNS,VofFluxEW,VofFluxNS,MoParExCo
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: FluxDiv,GradP
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: UFric,VFric,UWE,USN,VWE,VSN
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: matr
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: rhm,Uro,Vro
      INTEGER(kind=it4b):: num_iterations,Order2rd,Surface0
      REAL(KIND=dp):: final_res_norm
      REAL(KIND=dp):: Fe,Fw,Fn,Fs,ExDe,ExDw,ExDn,ExDs
      REAL(KIND=dp):: RoCeN,RoCeO,resi
      REAL(KIND=dp):: BetaP,BetaM,BetaW,Lamda,Fluxn0
  !   for particle tracking
      REAL(KIND=dp):: nupp,ropp,ug,vg,Reyp,Cd,mp,tp,Vrel
      allocate(MaFluxEW(ibeg:Isize+1,jbeg:Jsize,2))
      allocate(MaFluxNS(ibeg:Isize,jbeg:Jsize+1,2))
      allocate(VofFluxEW(ibeg:Isize+1,jbeg:Jsize,2))
      allocate(VofFluxNS(ibeg:Isize,jbeg:Jsize+1,2))
      allocate(CFluxEW(ibeg:Isize+1,jbeg:Jsize,2))
      allocate(CFluxNS(ibeg:Isize,jbeg:Jsize+1,2))
      allocate(ConFluxEW(ibeg:Isize+1,jbeg:Jsize,2))
      allocate(ConFluxNS(ibeg:Isize,jbeg:Jsize+1,2))
      allocate(DifFluxEW(ibeg:ibeg+Isize,jbeg:jbeg+Jsize-1,2))
      allocate(DifFluxNS(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize,2))
      allocate(EDFluxEW(ibeg:ibeg+Isize,jbeg:jbeg+Jsize-1,2))
      allocate(EDFluxNS(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize,2))
      allocate(ExEDFluxEW(ibeg:ibeg+Isize,jbeg:jbeg+Jsize-1,2))
      allocate(ExEDFluxNS(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize,2))

      allocate(FluxDiv(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1,2))
      allocate(MoParExCo(Isize,Jsize,2))
      allocate(MoParExIJU(TraPar%np,2))
      allocate(MoParExIJV(TraPar%np,2))
      allocate(GradP(Isize,Jsize,2))
      allocate(UWE(jbeg+Jsize-1,2))
      allocate(USN(ibeg+Isize-1,2))
      allocate(VWE(jbeg+Jsize-1,2))
      allocate(VSN(ibeg+Isize-1,2))
      allocate(matr(Isize,Jsize,5))
      allocate(rhm(Isize,Jsize))
      allocate(Uro(Isize,Jsize))
      allocate(Vro(Isize,Jsize))
      u => TVar%u
      v => TVar%v
      p => TVar%p
      Uvolf=>UCellO%vof
      Vvolf=>VCellO%vof
      do i=1,Isize
        do j=1,Jsize
          Uro(i,j)=(1.d0-Uvolf(i,j)-UCell%vofS(i,j))*roa/Roref+Uvolf(i,j)*      &
                                                                  row/Roref
          Vro(i,j)=(1.d0-Vvolf(i,j)-VCell%vofS(i,j))*roa/Roref+Vvolf(i,j)*      &
                                                                      row/Roref
        end do
      end do
      do i=1,TraPar%np
        call ParticlePosition(TraPar%Posp(i),UGrid,MoParExIJU(i,1),            &
                                                   MoParExIJU(i,2))
        call ParticlePosition(TraPar%Posp(i),VGrid,MoParExIJV(i,1),            &
                                                   MoParExIJV(i,2))
        if(MoParExIJU(i,1)==-1.or.MoParExIJU(i,2)==-1) then
          print*,itt
          print*,TraPar%Posp(i)%x,TraPar%Posp(i)%y
          print*,TraPar%dp(i)
          print*,
          print*,UGrid%x(1,1)-UGrid%dx(1,1)/2.d0
          print*,UGrid%x(1,1)+UGrid%dx(1,1)/2.d0
          print*,
          print*,TraPar%uvp(i)%u,TraPar%uvp(i)%v
          print*, 'problem with algorithm locating particles in UCell_PredictorUv_105'
          TraPar%uvp(i)%u=UParInlet
          MoParExIJU(i,1)=MoParExIJU(i-1,1)
          MoParExIJU(i,2)=MoParExIJU(i-1,2)
        end if
        if(MoParExIJV(i,1)==-1.or.MoParExIJV(i,2)==-1) then
          print*,itt
          print*,TraPar%Posp(i)%x,TraPar%Posp(i)%y
          print*,TraPar%dp(i)
          print*,VGrid%x(1,1)-VGrid%dx(1,1)/2.d0
          print*,VGrid%x(1,1)+VGrid%dx(1,1)/2.d0
          print*, 'problem with algorithm locating particles in VCell_PredictorUv_105'
          TraPar%uvp(i)%u=UParInlet
          MoParExIJV(i,1)=MoParExIJV(i-1,1)
          MoParExIJV(i,2)=MoParExIJV(i-1,2)
        end if
      end do
      BetaP=1.d0/(row/Roref)
      BetaM=1.d0/(roa/Roref)
!     Step 1: Calculate the convective coefficient
      if(itt==1) then
        call ModifiedConvectiveFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,    &
                                          CFluxEW,BoomCase%vs,1,0)
        call ModifiedConvectiveFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,    &
                                          CFluxNS,BoomCase%vs,0,1)
      else
 !       call SecondOrderConvectiveFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO, &
 !                                       CFluxEW,BoomCase%vs,1,0)
 !       call SecondOrderConvectiveFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO, &
 !                                       CFluxNS,BoomCase%vs,0,1)
        call HighOrderConvectiveFluxForXDir(PGrid,UGrid,VGrid,UCellO,VCellO,   &
                                       BoomCase%vs,CFluxEW)
        call HighOrderConvectiveFluxForYDir(PGrid,UGrid,VGrid,UCellO,VCellO,     &
                                       BoomCase%vs,CFluxNS)
      end if
!      call HighOrderDensityBasedConvFluxXDir(PGrid,UGrid,VGrid,PCell,UCell,    &
!                                      VCell,Uro,Vro,ConfluxEW,MaFluxEW,dt)
!      call HighOrderDensityBasedConvFluxYDir(PGrid,UGrid,VGrid,PCell,UCell,    &
!                                      VCell,Uro,Vro,ConfluxNS,MaFluxNS,dt)
      call FaceDensityFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,MaFluxEW,    &
                                             VofFluxEW,BoomCase%vs,dt,1,0)
      call FaceDensityFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,MaFluxNS,    &
                                             VofFluxNS,BoomCase%vs,dt,0,1)
      call DensityBasedConvectiveFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,  &
                                      ConFluxEW,VofFluxEW,BoomCase%vs,1,0)
      call DensityBasedConvectiveFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,  &
                                      ConFluxNS,VofFluxNS,BoomCase%vs,0,1)
!     Step 2: Calculate the diffusive coefficient
      call DiffusiveFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,DifFluxEW,     &
                                           EDFluxEW,ExEDFluxEW,BoomCase%vs,1,0)
      call DiffusiveFlux(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,DifFluxNS,     &
                                           EDFluxNS,ExEDFluxNS,BoomCase%vs,0,1)
!     Step 3: Calculate source term coefficient as wall function
!     for U Cell
      do i = 1,Isize-1
        do j = 1,Jsize
          GradP(i,j,1)=0.d0
          FluxDiv(i,j,1)=0.d0
          if(UCell%Posnu(i,j)/=-1) then
            Fe=ConFluxEW(i+1,j,1)
            Fw=ConFluxEW(i,j,1)
            Fn=ConFluxNS(i,j+1,1)
            Fs=ConFluxNS(i,j,1)
            RoCeO=(1.d0-Uvolf(i,j)-UCell%vofS(i,j))*roa/Roref+Uvolf(i,j)*      &
                                                                  row/Roref
            RoCeN=RoCeO-dt*(MaFluxEW(i+1,j,1)-MaFluxEW(i,j,1)+                 &
                            MaFluxNS(i,j+1,1)-MaFluxNS(i,j,1))/                &
                                              (UGrid%dx(i,j)*UGrid%dy(i,j))
            if(RoCeN>(1.d0+epsi-UCell%vofS(i,j))*roa/Roref.and.                &
               RoCeN<(1.d0-epsi-UCell%vofS(i,j))*row/Roref.and.                &
               RoCeO>(1.d0+epsi-UCell%vofS(i,j))*roa/Roref.and.                &
               RoCeO<(1.d0-epsi-UCell%vofS(i,j))*row/Roref.and.                &
               UCell%vofS(i,j)<epsi)then
              Fluxn0=(Fe-Fw+Fn-Fs)/RoCeN+u(i,j)*(1.d0-RoCeO/RoCeN)/dt*         &
                                                UGrid%dx(i,j)*UGrid%dy(i,j)
            else
              Fluxn0=(CFluxEW(i+1,j,1)-CFluxEW(i,j,1)+CFluxNS(i,j+1,1)-        &
                      CFluxNS(i,j,1))/(1.d0-UCell%vofS(i,j))

            end if
            if(((RoCeN>(1.d0+epsi-UCell%vofS(i,j))*roa/Roref.and.              &
                RoCeN<(1.d0-epsi-UCell%vofS(i,j))*row/Roref.and.               &
                RoCeO>(1.d0+epsi-UCell%vofS(i,j))*roa/Roref.and.               &
                RoCeO<(1.d0-epsi-UCell%vofS(i,j))*row/Roref.and.               &
                UCell%vofS(i,j)<epsi).or.                                      &
                RoCeN<(1.d0+epsi-UCell%vofS(i,j))*roa/Roref.or.                &
                RoCeN>(1.d0-epsi-UCell%vofS(i,j))*row/Roref))  then
              Order2rd=1
            else
              Order2rd=0
            end if
            Order2rd=1
            if(itt>1.and.Order2rd==1) then
              FluxDiv(i,j,1)=1.5d0*Fluxn0-0.5d0*Flux_n1(i,j,1)
            else
              FluxDiv(i,j,1)=Fluxn0
            end if
          ! for Cell with all open cell faces equaled to 0
          ! set up fluxdiv(i,j,1)=0.
            Surface0=0
            if(UCell%WEdge_Area(i,j)<epsi) Surface0=Surface0+1
            if(UCell%EEdge_Area(i,j)<epsi) Surface0=Surface0+1
            if(UCell%SEdge_Area(i,j)<epsi) Surface0=Surface0+1
            if(UCell%NEdge_Area(i,j)<epsi) Surface0=Surface0+1
            if(Surface0>=3) then
              Surface0=0
              FluxDiv(i,j,1)=0.d0
            end if
            FluxDiv(i,j,1)=FluxDiv(i,j,1)-u(i,j)*(1.d0-UCellO%vofS(i,j))/      &
                         (1.d0-UCell%vofS(i,j))*UGrid%dx(i,j)*UGrid%dy(i,j)/dt
            Flux_n1(i,j,1)=Fluxn0
          end if
        end do
      end do
      call Mixing_Procedure(UGrid,UCell,1,0,FluxDiv(:,:,1))
      do i=1,Isize-1
        do j=1,Jsize
          if(UCell%Posnu(i,j)/=-1) then
            ExDe = ExEDFluxEW(i+1,j,1)
            ExDw = ExEDFluxEW(i,j,1)
            ExDn = ExEDFluxNS(i,j+1,1)
            ExDs = ExEDFluxNS(i,j,1)
            FluxDiv(i,j,1)=FluxDiv(i,j,1)-                                     &
                          BetaVis*((DifFluxEW(i+1,j,1)*(u(i+1,j)-u(i,j))-      &
                                    DifFluxEW(i,j,1)*(u(i,j)-u(i-1,j)))-       &
                                   (DifFluxNS(i,j+1,1)*(u(i,j+1)-u(i,j))-      &
                                    DifFluxNS(i,j,1)*(u(i,j)-u(i,j-1))))/      &
                                   (1.d0-UCell%vofS(i,j))+                     &
                                   (-ExDe+ExDw-ExDn+ExDs)/(1.d0-UCell%vofS(i,j))
            if(isnan(FluxDiv(i,j,1)).or.dabs(Fluxdiv(i,j,1))>1.d10) then
              print*,'Problem with U-Velocity'
              print*,itt
              print*,i,j
              print*,Pred%u(i,j),UCell%vofS(i,j)
              print*,'flux'
              pause'PredictorUV_156'
            end if
          end if
        end do
      end do
 !    Adding pressure gradient
      do i = ibeg,Isize+ibeg-2
        do j = jbeg,Jsize+jbeg-1
          if(UCell%Posnu(i,j)/=-1) then
      !   Add pressure gradient
            if((PCellO%vof(i,j)>=0.5d0.and.PCellO%vofS(i,j)<epsi).or.          &
               (PCellO%phi(i,j)<0.d0.and.PCellO%vofS(i,j)>=epsi)) then
              if((PCellO%vof(i+1,j)<0.5d0.and.PCellO%vofS(i+1,j)<epsi).or.     &
                 (PCellO%phi(i+1,j)>=0.d0.and.PCellO%vofS(i+1,j)>=epsi))then
                Lamda=dabs(PCellO%phi(i,j))/(dabs(PCellO%phi(i,j))+            &
                                             dabs(PCellO%phi(i+1,j))+tol)
                BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP
                GradP(i,j,1)=UGrid%dy(i,j)*(p(i+1,j)-p(i,j))*BetaM*BetaP/BetaW
              else
                BetaW = BetaP
                GradP(i,j,1)=UGrid%dy(i,j)*(p(i+1,j)-p(i,j))*BetaW
              end if
            else
              if((PCellO%vof(i+1,j)>=0.5d0.and.PCellO%vofS(i+1,j)<epsi).or.    &
                 (PCellO%phi(i+1,j)<0.d0.and.PCellO%vofS(i+1,j)>=epsi))then
                Lamda=dabs(PCellO%phi(i,j))/(dabs(PCellO%phi(i,j))+            &
                                             dabs(PCellO%phi(i+1,j))+tol)
                BetaW=Lamda*BetaP+(1.d0-Lamda)*BetaM
                GradP(i,j,1)=UGrid%dy(i,j)*(p(i+1,j)-p(i,j))*BetaM*BetaP/BetaW
              else
                BetaW=BetaM
                GradP(i,j,1)=UGrid%dy(i,j)*(p(i+1,j)-p(i,j))*BetaW
              end if
            end if
            Gradp(i,j,1)=0.d0
            FluxDiv(i,j,1)=FluxDiv(i,j,1)!+GradP(i,j,1)
          end if
        end do
      end do

 !    Add the particle momentum exchange to flux difference for UCell
      do i=1,TraPar%np
        ii=MoParExIJU(i,1)
        jj=MoParExIJU(i,2)
        nupp=nuw*UCell%vof(ii,jj)/(1.d0-UCell%vofS(ii,jj))+                    &
          nua*(1.d0-UCell%vof(ii,jj)-UCell%vofS(ii,jj))/(1.d0-UCell%vofS(ii,jj))
        ropp=row*UCell%vof(ii,jj)/(1.d0-UCell%vofS(ii,jj))+                    &
          roa*(1.d0-UCell%vof(ii,jj)-UCell%vofS(ii,jj))/(1.d0-UCell%vofS(ii,jj))
        mp=Rop*1.d0/6.d0*pi*TraPar%dp(i)**3.d0
        ug=TVar%u(ii,jj)
        vg=TVar%v(ii,jj)
        VRel=dsqrt((TraPar%uvp(i)%u-ug)**2.d0+(TraPar%uvp(i)%v-vg)**2.d0)
        if(VRel<1.d-7) VRel=1.d-7
        Reyp=TraPar%dp(i)*VRel/nupp
        Cd=Drag(Reyp)
        tp=4.d0/3.d0*TraPar%dp(i)*rop/ropp/Cd
        if(ii/=-1.and.jj/=-1) then
          fluxDiv(ii,jj,1)=fluxDiv(ii,jj,1)+mp/tp*TraPar%uvp(i)%u*VRel/        &
                          (TVar%Roref*TVar%Uref**2.d0/UGrid%Lref)/zp/Uro(ii,jj)
          MoParExCo(ii,jj,1)=MoParExCo(ii,jj,1)+mp/tp*VRel/                    &
                          (TVar%Roref*TVar%Uref/UGrid%Lref)/zp/Uro(ii,jj)
          if(dabs(MoParExco(ii,jj,1))>1.d10.or.isnan(MoParExco(ii,jj,1)).or.   &
             isnan(FluxDiv(ii,jj,1))) then
            print*,itt
            print*,ii,jj
            print*,mp,tp
            print*,VRel,(TVar%Roref*TVar%Uref/UGrid%Lref)/zp/Uro(ii,jj)
            pause'Momentum particle 280'
          end if
        end if
      end do
     !   For V-Cell
      do i=1,Isize
        do j=1,Jsize-1
          FluxDiv(i,j,2)=0.d0
          if(VCell%Posnu(i,j)/=-1) then
            Fe=ConFluxEW(i+1,j,2)
            Fw=ConFluxEW(i,j,2)
            Fn=ConFluxNS(i,j+1,2)
            Fs=ConFluxNS(i,j,2)
            RoCeO=(1.d0-Vvolf(i,j)-VCell%vofS(i,j))*roa/Roref+Vvolf(i,j)*      &
                                                                      row/Roref
            RoCeN=RoCeO-dt*(MaFluxEW(i+1,j,2)-MaFluxEW(i,j,2)+                 &
                            MaFluxNS(i,j+1,2)-MaFluxNS(i,j,2))/                &
                                              (VGrid%dx(i,j)*VGrid%dy(i,j))
            if(RoCeN>(1.d0+epsi-VCell%vofS(i,j))*roa/Roref.and.                &
               RoCeN<(1.d0-epsi-VCell%vofS(i,j))*row/Roref.and.                &
               RoCeO>(1.d0+epsi-VCell%vofS(i,j))*roa/Roref.and.                &
               RoCeO<(1.d0-epsi-VCell%vofS(i,j))*row/Roref.and.                &
               VCell%vofS(i,j)<epsi) then
              Fluxn0=(Fe-Fw+Fn-Fs)/RoCeN+v(i,j)*(1.d0-RoCeO/RoCeN)/dt*         &
                                            VGrid%dx(i,j)*VGrid%dy(i,j)
            else
              Fluxn0=(CFluxEW(i+1,j,2)-CFluxEW(i,j,2)+CFluxNS(i,j+1,2)-        &
                                       CFluxNS(i,j,2))/(1.d0-VCell%VofS(i,j))
            end if
            if(((RoCeN>(1.d0+epsi-VCell%vofS(i,j))*roa/Roref.and.              &
                RoCeN<(1.d0-epsi-VCell%vofS(i,j))*row/Roref.and.               &
                RoCeO>(1.d0+epsi-VCell%vofS(i,j))*roa/Roref.and.               &
                RoCeO<(1.d0-epsi-VCell%vofS(i,j))*row/Roref.and.               &
                VCell%vofS(i,j)<epsi).or.                                      &
                RoCeN<(1.d0+epsi-VCell%vofS(i,j))*roa/Roref.or.                &
                RoCeN>(1.d0-epsi-VCell%vofS(i,j))*row/Roref)) then
              Order2rd=1
            else
              Order2rd=0
            end if
            Order2rd=1
      !     V Cell
            if(itt>1.and.Order2rd==1) then
              FluxDiv(i,j,2)=1.5d0*Fluxn0-0.5d0*Flux_n1(i,j,2)
            else
              FluxDiv(i,j,2)=Fluxn0
            end if
            Surface0=0
            if(VCell%WEdge_Area(i,j)<epsi) Surface0=Surface0+1
            if(VCell%EEdge_Area(i,j)<epsi) Surface0=Surface0+1
            if(VCell%SEdge_Area(i,j)<epsi) Surface0=Surface0+1
            if(VCell%NEdge_Area(i,j)<epsi) Surface0=Surface0+1
            if(Surface0>=3) then
              Surface0=0
              FluxDiv(i,j,2)=0.d0
            end if
            FluxDiv(i,j,2)=FluxDiv(i,j,2)-v(i,j)*(1.d0-VCellO%VofS(i,j))/      &
                       (1.d0-VCell%VofS(i,j))*VGrid%dx(i,j)*VGrid%dy(i,j)/dt
            Flux_n1(i,j,2)=Fluxn0
           ! if(itt>195.and.(i==306.and.j==116)) then
           !   print*,'test flux'
           !   print*,itt
           !   print*,FluxDiv(i,j,2)
           !   print*,v(i,j)
           !   print*,v(i,j)*(1.d0-VCellO%VofS(i,j))/                           &
           !            (1.d0-VCell%VofS(i,j))*VGrid%dx(i,j)*VGrid%dy(i,j)/dt
           !   print*,
           ! end if
          end if
        end do
      end do
      call Mixing_Procedure(VGrid,VCell,0,1,FluxDiv(:,:,2))
      do i=1,Isize
        do j=1,Jsize-1
          if(VCell%Posnu(i,j)/=-1) then
            ExDe = ExEDFluxEW(i+1,j,2)
            ExDw = ExEDFluxEW(i,j,2)
            ExDn = ExEDFluxNS(i,j+1,2)
            ExDs = ExEDFluxNS(i,j,2)
            FluxDiv(i,j,2)=FluxDiv(i,j,2)-                                     &
                  BetaVis*((DifFluxEW(i+1,j,2)*(v(i+1,j)-v(i,j))-              &
                            DifFluxEW(i,j,2)*(v(i,j)-v(i-1,j)))-               &
                           (DifFluxNS(i,j+1,2)*(v(i,j+1)-v(i,j))-              &
                            DifFluxNS(i,j,2)*(v(i,j)-v(i,j-1))))/              &
                           (1.d0-VCell%vofS(i,j))+                             &
                           (-ExDe+ExDw-ExDn+ExDs)/(1.d0-VCell%vofS(i,j))
            if(VCell%VofS(i,j)>epsi.and.VCell%VofS(i,j)<1.d0-epsi) then
              !*dabs(VCell%nx(i,j))
              FluxDiv(i,j,2)=FluxDiv(i,j,2)-((1.d0-VCell%vof(i,j)/             &
                            (1.d0-VCell%vofS(i,j)))*nua/nuref+VCell%vof(i,j)/  &
                            (1.d0-VCell%vofS(i,j))*nuw/nuref)/Rey*             &
                            VCell%WlLh(i,j)/VCell%delh(i,j)*BoomCase%vs
            end if
            if(isnan(FluxDiv(i,j,2)).or.dabs(Fluxdiv(i,j,2))>1.d20) then
              print*,'Problem with V-Velocity'
              print*,itt
              print*,i,j
              print*,VCell%vofS(i,j)
              print*,VCell%Posnu(i,j)
              print*,
              print*,FluxDiv(i,j,2)
              pause'PredictorUV_218'
            end if
          end if
        end do
      end do
      do i=1,Isize
        do j=1,Jsize-1
          if(VCell%Posnu(i,j)/=-1) then
      !     Add pressure gradient
            if((PCellO%vof(i,j)>=0.5d0.and.PCellO%vofS(i,j)<epsi).or.          &
                (PCellO%phi(i,j)<0.d0.and.PCellO%vofS(i,j)>=epsi)) then
              if((PCellO%vof(i,j+1)<0.5d0.and.PCellO%vofS(i,j+1)<epsi).or.     &
                 (PCellO%phi(i,j+1)>=0.d0.and.PCellO%vofS(i,j+1)>=epsi)) then
                Lamda=dabs(PCellO%phi(i,j))/(dabs(PCellO%phi(i,j))+            &
                                             dabs(PCellO%phi(i,j+1))+tol)
                BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP
                GradP(i,j,2)=VGrid%dx(i,j)*(p(i,j+1)-p(i,j))*BetaP*BetaM/BetaW
              else
                BetaW=BetaP
                GradP(i,j,2)=VGrid%dx(i,j)*(p(i,j+1)-p(i,j))*BetaW
              end if
            else
              if((PCellO%vof(i,j+1)>=0.5d0.and.PCellO%vofS(i,j+1)<epsi).or.    &
                 (PCellO%phi(i,j+1)<0.d0.and.PCellO%vofS(i,j+1)>=epsi))  then
                Lamda=dabs(PCellO%phi(i,j))/(dabs(PCellO%phi(i,j))+            &
                                             dabs(PCellO%phi(i,j+1))+tol)
                BetaW=Lamda*BetaP+(1.d0-Lamda)*BetaM
                GradP(i,j,2)=VGrid%dx(i,j)*(p(i,j+1)-p(i,j))*BetaP*BetaM/BetaW
              else
                BetaW=BetaM
                GradP(i,j,2)=VGrid%dx(i,j)*(p(i,j+1)-p(i,j))*BetaW
              end if
            end if
            GradP(i,j,2)=0.d0
            FluxDiv(i,j,2)=FluxDiv(i,j,2)!+GradP(i,j,2)            ! Contribution of pressure gradient
          end if
        end do
      end do
    ! Add the particle momentum exchange to flux difference for UCell
      do i=1,TraPar%np
        ii=MoParExIJU(i,1)
        jj=MoParExIJU(i,2)
        nupp=nuw*VCell%vof(ii,jj)/(1.d0-VCell%vofS(ii,jj))+                    &
          nua*(1.d0-VCell%vof(ii,jj)-VCell%vofS(ii,jj))/(1.d0-VCell%vofS(ii,jj))
        ropp=row*VCell%vof(ii,jj)/(1.d0-VCell%vofS(ii,jj))+                    &
          roa*(1.d0-VCell%vof(ii,jj)-VCell%vofS(ii,jj))/(1.d0-VCell%vofS(ii,jj))
        mp=Rop*1.d0/6.d0*pi*TraPar%dp(i)**3.d0
        ug=TVar%u(ii,jj)
        vg=TVar%v(ii,jj)
        VRel=dsqrt((TraPar%uvp(i)%u-ug)**2.d0+(TraPar%uvp(i)%v-vg)**2.d0)
        if(VRel<1.d-7) VRel=1.d-7
        Reyp=TraPar%dp(i)*VRel/nupp
        Cd=Drag(Reyp)
        tp=4.d0/3.d0*TraPar%dp(i)*rop/ropp/Cd
        if(ii/=-1.and.jj/=-1) then
          FluxDiv(ii,jj,2)=FluxDiv(ii,jj,2)+mp/tp*TraPar%uvp(i)%v*VRel/        &
                          (TVar%Roref*TVar%Uref**2.d0/VGrid%Lref)/zp/Vro(ii,jj)
          MoParExCo(ii,jj,2)=MoParExCo(ii,jj,2)+mp/tp*VRel/(TVar%Roref*        &
                           TVar%Uref/VGrid%Lref)/zp/Vro(ii,jj)
          if(dabs(MoParExco(ii,jj,2))>1.d10.or.isnan(fluxDiv(ii,jj,2)).or.     &
             isnan(MoParExco(ii,jj,2))) then
            print*,fluxDiv(ii,jj,2)
            print*,TraPar%mp(i),TraPar%tp(i)
            print*,MoParExco(ii,jj,2)
            print*,TraPar%mp(i)/TraPar%tp(i)*TraPar%uvp(i)%v*TraPar%VRelG(i)/  &
                          (TVar%Roref*TVar%Uref**2.d0/VGrid%Lref)/zp/Vro(ii,jj)
            print*,TraPar%mp(i)/TraPar%tp(i)/TraPar%VRelG(i)/                  &
                  (TVar%Roref*TVar%Uref/VGrid%Lref)/zp/Vro(ii,jj)
            print*,
            print*,TraPar%mp(i),TraPar%tp(i),TraPar%VRelG(i)
            print*,TVar%Roref*TVar%Uref/VGrid%Lref,zp,Vro(ii,jj)
            pause'Momentum Particle 425'
          end if
        end if
      end do
    ! Solving for UCell
      nullify(Uvolf,Vvolf)
      Uvolf=>UCell%vof
      Vvolf=>VCell%vof
      call DiffusiveFlux(PGrid,UGrid,VGrid,PCell,UCell,VCell,DifFluxEW,        &
                                           EDFluxEW,ExEDfluxEW,BoomCase%vs,1,0)
      call DiffusiveFlux(PGrid,UGrid,VGrid,PCell,UCell,VCell,DifFluxNS,        &
                                           EDFluxNS,ExEDfluxEW,BoomCase%vs,0,1)
      nullify(Uvolf,Vvolf)
      call SetBasicSolver(solver,precond)
      call SetMatrix(A,parcsr_A,UGrid,UCell,DifFluxEW,DifFluxNS,EDFluxEW,      &
                                EDFluxNS,PU,UWE,USN,matr,dt,1,0)
      call SetVectors(b,x,par_b,par_x,PGrid,UGrid,UCell,UWE,USN,FluxDiv(:,:,1),&
                                                     rhm,dt,1,0)
      call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
      call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    ! Run info - needed logging turned on
      call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      call DeltaGetValues(x,UCell,Pred%u,1,0)
      do i = 1,Isize-1
        do j = 1,Jsize
          Pred%u(i,j)=Pred%u(i,j)+dt*GradP(i,j,1)/UGrid%dx(i,j)/UGrid%dy(i,j)
        end do
      end do

      call HYPRE_IJMatrixDestroy(A,ierr)
      call HYPRE_IJVectorDestroy(b,ierr)
      call HYPRE_IJVectorDestroy(x,ierr)
      call HYPRE_BoomerAMGDestroy(precond,ierr)
      call HYPRE_ParCSRPCGDestroy(solver,ierr)
    ! For VCell
      call SetBasicSolver(solver,precond)
      call SetMatrix(A,parcsr_A,VGrid,VCell,DifFluxEW,DifFluxNS,EDFluxEW,      &
                                           EDFluxNS,PV,VWE,VSN,matr,dt,0,1)
      call SetVectors(b,x,par_b,par_x,PGrid,VGrid,VCell,VWE,VSN,FluxDiv(:,:,2),&
                                                                   rhm,dt,0,1)
      call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
      call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
    ! Run info - needed logging turned on
      call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      call DeltaGetValues(x,VCell,Pred%v,0,1)
      do i = 1,Isize
        do j = 1,Jsize-1
          Pred%v(i,j)=Pred%v(i,j)+dt*GradP(i,j,2)/VGrid%dx(i,j)/VGrid%dy(i,j)
        end do
      end do

      call HYPRE_IJMatrixDestroy(A,ierr)
      call HYPRE_IJVectorDestroy(b,ierr)
      call HYPRE_IJVectorDestroy(x,ierr)
      call HYPRE_BoomerAMGDestroy(precond,ierr)
      call HYPRE_ParCSRPCGDestroy(solver,ierr)
      PU%Dp(ibeg-ight,:) = PU%Dp(ibeg,:)
      PU%Dp(ibeg+Isize-1,:) = PU%Dp(ibeg+Isize-2,:)
      PU%Dp(ibeg+Isize-1+ight,:) = PU%Dp(ibeg+Isize-1,:)
      PU%Dp(:,jbeg-jght) = PU%Dp(:,jbeg)
      PU%Dp(:,jbeg+Jsize-1+jght) = PU%Dp(:,jbeg+Jsize-1)

      PV%Dp(ibeg-ight,:) = PV%Dp(ibeg,:)
      PV%Dp(ibeg+Isize-1+ight,:) = PV%Dp(ibeg+Isize-1,:)
      PV%Dp(:,jbeg-jght) = PV%Dp(:,jbeg)
      PV%Dp(:,jbeg+Jsize-1) = PV%Dp(:,jbeg+Jsize-2)
      PV%Dp(:,jbeg+Jsize-1+jght) = PV%Dp(:,jbeg+Jsize-1)
      call PredictorVelocityBoundaryCondition(Pred,TVar)
      deallocate(Uro,Vro)
      deallocate(MaFluxEW)
      deallocate(MaFluxNS)
      deallocate(VofFluxEW)
      deallocate(VofFluxNS)
      deallocate(CFluxEW)
      deallocate(CFluxNS)
      deallocate(ConFluxEW)
      deallocate(ConFluxNS)
      deallocate(DifFluxEW,EDFluxEW,ExEDFluxEW)
      deallocate(DifFluxNS,EDFluxNS,ExEDFluxNS)
      deallocate(UWE)
      deallocate(USN)
      deallocate(VWE)
      deallocate(VSN)
      deallocate(FluxDiv)
      deallocate(GradP)
      deallocate(matr,rhm)
      nullify(u)
      nullify(v)
      nullify(p)
    End subroutine Predictor_UV

!   Calculate normal flux for both normal cell and cut cell
    Subroutine SetBasicSolver(solver,precond)
        IMPLICIT NONE
        INTEGER*8,INTENT(INOUT):: solver
        INTEGER*8,INTENT(INOUT),optional:: precond
!       Set up and use a solver
        Call HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD,solver,ierr)
!       Set some PARAMETERs
        Call HYPRE_ParCSRPCGSetMaxIter(solver,50,ierr)
        Call HYPRE_ParCSRPCGSetTol(solver,1.0d-20,ierr)
        Call HYPRE_ParCSRPCGSetTwoNorm(solver,0,ierr)
!        Call HYPRE_ParCSRPCGSetPrintLevel(solver,2,ierr)
        Call HYPRE_ParCSRPCGSetLogging(solver,1,ierr)
!       Now set up the AMG preconditioner and specify any PARAMETERs
        If(present(precond)) then
          Call HYPRE_BoomerAMGCreate(precond,ierr)
!        Set some PARAMETERs
!        Print less solver info since a preconditioner
!          Call HYPRE_BoomerAMGSetPrintLevel(precond,1,ierr);
!        Falgout coarsening
          Call HYPRE_BoomerAMGSetCoarsenTYPE(precond,6,ierr)
!        SYMMETRIC G-S/Jacobi hybrid relaxation
          Call HYPRE_BoomerAMGSetRelaxTYPE(precond,6,ierr)
!        Sweeeps on each level
          Call HYPRE_BoomerAMGSetNumSweeps(precond,1,ierr)
!        conv. tolerance
          Call HYPRE_BoomerAMGSetTol(precond,0.0d0,ierr)
!        do only one iteration!
          Call HYPRE_BoomerAMGSetMaxIter(precond,1,ierr)
!        set amg as the pcg preconditioner
!         precond_id = 2
          Call HYPRE_ParCSRPCGSetPrecond(solver,2,precond,ierr)
        End if
    End subroutine SetBasicSolver

    subroutine SetMatrix(A,parcsr_A,TGrid,TCell,DifFluxEW,DifFluxNS,EDFluxEW,  &
                                                   EDFluxNS,PUV,CWE,CSN,matr,dt,iu,iv)
        IMPLICIT NONE
        INTEGER*8,INTENT(INOUT):: A,parcsr_A
        TYPE(Grid),INTENT(IN):: TGrid
        TYPE(Cell),INTENT(IN):: TCell
        REAL(KIND=dp),INTENT(IN):: dt
        REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(IN):: DifFluxEW,     &
                                                   DifFluxNS,EDFluxEW,EDFluxNS
        TYPE(PoissonCoefficient),INTENT(INOUT):: PUV
        REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(INOUT):: CWE,CSN
        REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT):: matr
        INTEGER,INTENT(IN):: iu,iv
        INTEGER(kind=it4b):: nnz,ictr,ilower,iupper,cols(0:4)
        INTEGER(kind=it4b):: i,j
        REAL(KIND=dp):: aP,aE,aW,aN,aS,De,Dw,Dn,Ds,Sp,VofAF
        REAL(KIND=dp):: Fe,Fw,Fn,Fs,Fep,Fem,Fwp,Fwm,Fnp,Fnm,Fsp,Fsm
        REAL(KIND=dp):: values(0:4)
        ilower = 0
        iupper = TCell%ExtCell
      ! Create and Set up matrix
        call HYPRE_IJMatrixCreate(MPI_COMM_WORLD,ilower,iupper,ilower,iupper,  &
                                                                        A,ierr)
        call HYPRE_IJMatrixSetObjectTYPE(A,HYPRE_PARCSR,ierr)
        call HYPRE_IJMatrixInitialize(A,ierr)
        do i = ibeg,Isize-iu
          do j = jbeg,Jsize-iv
            if(TCell%Posnu(i,j)/=-1) then
              matr(i,j,:)=0.d0
              VofAF=1.d0-TCell%vofS(i,j)
              De=0.d0;Dw=0.d0;Dn=0.d0;Ds=0.d0
              if(iu==1) then ! for UCell
                De=(1.d0-BetaVis)*DifFluxEW(i+1,j,1)/VofAF
                Dw=(1.d0-BetaVis)*DifFluxEW(i,j,1)/VofAF
                Dn=(1.d0-BetaVis)*DifFluxNS(i,j+1,1)/VofAF
                Ds=(1.d0-BetaVis)*DifFluxNS(i,j,1)/VofAF
                Fep=EDfluxEW(i+1,j,1)*(1.d0-TCell%EtaE(i,j))/VofAF
                Fem=EDfluxEW(i+1,j,1)*TCell%EtaE(i,j)/VofAF
                Fwp=EDfluxEW(i,j,1)*TCell%EtaE(i-1,j)/VofAF
                Fwm=EDfluxEW(i,j,1)*(1.d0-TCell%EtaE(i-1,j))/VofAF
                Fnp=EDfluxNS(i,j+1,1)*(1.d0-TCell%EtaN(i,j))/VofAF
                Fnm=EDfluxNS(i,j+1,1)*TCell%EtaN(i,j)/VofAF
                Fsp=EDfluxNS(i,j,1)*TCell%EtaN(i,j-1)/VofAF
                Fsm=EDfluxNS(i,j,1)*(1.d0-TCell%EtaN(i,j-1))/VofAF
                if(j==Jsize)Dn = 0.d0
                if(i==Isize-iu) De = 0.d0
                if(i==ibeg)Dw = 0.d0
              elseif(iv==1) then ! for VCell
                De=(1.d0-BetaVis)*DifFluxEW(i+1,j,2)/VofAF
                Dw=(1.d0-BetaVis)*DifFluxEW(i,j,2)/VofAF
                Dn=(1.d0-BetaVis)*DifFluxNS(i,j+1,2)/VofAF
                Ds=(1.d0-BetaVis)*DifFluxNS(i,j,2)/VofAF
                Fep=EDfluxEW(i+1,j,2)*(1.d0-TCell%EtaE(i,j))/VofAF
                Fem=EDfluxEW(i+1,j,2)*TCell%EtaE(i,j)/VofAF
                Fwp=EDfluxEW(i,j,2)*TCell%EtaE(i-1,j)/VofAF
                Fwm=EDfluxEW(i,j,2)*(1.d0-TCell%EtaE(i-1,j))/VofAF
                Fnp=EDfluxNS(i,j+1,2)*(1.d0-TCell%EtaN(i,j))/VofAF
                Fnm=EDfluxNS(i,j+1,2)*TCell%EtaN(i,j)/VofAF
                Fsp=EDfluxNS(i,j,2)*TCell%EtaN(i,j-1)/VofAF
                Fsm=EDfluxNS(i,j,2)*(1.d0-TCell%EtaN(i,j-1))/VofAF
                if(j==Jsize-iv) Dn=0.d0
              end if
              aE=De;aW=Dw;aN=Dn;aS=Ds
              if(i==ibeg)CWE(j,1)=aW
              if(i==Isize-iu)CWE(j,2)=aE
              if(j==jbeg)CSN(i,1)=aS
              if(j==Jsize-iv)CSN(i,2)=aN
              aP=TGrid%dx(i,j)*TGrid%dy(i,j)/dt+aE+aW+aN+aS+                    &
                 ((1.d0-TCell%vof(i,j)/(1.d0-TCell%vofS(i,j)))*nua/nuref+       &
                   TCell%vof(i,j)/(1.d0-TCell%vofS(i,j))*nuw/nuref)/Rey*        &
                   TCell%WlLh(i,j)/TCell%delh(i,j)/VofAF
              aP=aP+Fep-Fwp+Fnp-Fsp
              aE=aE-Fem
              aW=aW+Fwm
              aN=aN-Fnm
              aS=aS+Fsm
          !   PUV%Dp(i,j)=dt !1.d0/(aP-aE-aW-aN-aS)
              PUV%Dp(i,j)=TGrid%dx(i,j)*TGrid%dy(i,j)/(aP-aE-aW-aN-aS)
              if(i==3001.and.j==125) then
                print*,'Test matrix coefficients'
                print*,aP
                print*,((1.d0-TCell%vof(i,j)/(1.d0-TCell%vofS(i,j)))*nua/nuref+       &
                   TCell%vof(i,j)/(1.d0-TCell%vofS(i,j))*nuw/nuref)/Rey*        &
                   TCell%WlLh(i,j)/TCell%delh(i,j)/VofAF
                print*,'End test matrix'
              end if
              ictr=TCell%Posnu(i,j)
              nnz=0
              values=0.d0
              cols=0
          !   Bottom of current cell
              if(j>jbeg) then
                if(TCell%Posnu(i,j-1)/=-1) then
                  cols(nnz)=TCell%Posnu(i,j-1)
                  values(nnz)=-aS
                  matr(i,j,1)=values(nnz)
                  nnz=nnz+1
                end if
              end if
          !   West of current cell
              if(i>ibeg) then
                if(TCell%Posnu(i-1,j)/=-1) then
                  cols(nnz)=TCell%Posnu(i-1,j)
                  values(nnz)=-aW
                  matr(i,j,2)=values(nnz)
                  nnz=nnz+1
                end if
              end if
          ! Set the diagonal cell
              cols(nnz)=TCell%Posnu(i,j)
              values(nnz)=aP
              if(isnan(values(nnz)).or.dabs(values(nnz))>1.d20) then
                print*,i,j,values(nnz)
                print*,iu
                print*,
                print*,aE,aW
                print*,aN,aS
                print*,
                print*,Dn,Fnm
                print*,DifFluxNS(i,j+1,1)
                print*,'test UVCell'
                print*,TCell%SyN(i,j),TCell%SyN(i,j+1),TCell%SyN(i,j-1)
                print*,TCell%Cell_Cent(i,j,2),TCell%Cell_Cent(i,j+1,2),TCell%Cell_Cent(i,j-1,2)
                print*,
                print*,Sp
                pause 'predictoruv_760'
              end if
              matr(i,j,3)=values(nnz)
              nnz = nnz+1
          ! East of current cell
              if(i<Isize) then
                if(TCell%Posnu(i+1,j)/=-1) then
                  cols(nnz) = TCell%Posnu(i+1,j)
                  values(nnz) = -aE
                  matr(i,j,4)=values(nnz)
                  nnz = nnz+1
                end if
              end if
          ! North of current cell
              if(j<Jsize) then
                if(TCell%Posnu(i,j+1)/=-1) then
                  cols(nnz) = TCell%Posnu(i,j+1)
                  values(nnz) = -aN
                  matr(i,j,5)=values(nnz)
                  nnz = nnz+1
                end if
              end if
              call HYPRE_IJMatrixSetValues(A,1,nnz,ictr,cols,values,ierr)
            end if
          end do
        end do
        call HYPRE_IJMatrixAssemble(A,ierr)
        call HYPRE_IJMatrixGetObject(A,parcsr_A,ierr)
    end subroutine SetMatrix

    subroutine SetVectors(b,x,par_b,par_x,PGrid,TGrid,TCell,CWE,CSN,IJFlux,    &
                                                            rhm,dt,iu,iv)
        INTEGER*8:: b,x,par_b,par_x
        INTEGER,INTENT(IN):: iu,iv
        TYPE(Grid),INTENT(IN):: PGrid,TGrid
        TYPE(Cell),INTENT(IN):: TCell
        REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(IN):: CWE,CSN
        REAL(KIND=dp),DIMENSION(:,:),INTENT(IN):: IJFlux
        REAL(KIND=dp),DIMENSION(:,:),INTENT(INOUT):: rhm
        real(dp),INTENT(IN):: dt
        INTEGER(kind=it4b):: i,j
        INTEGER:: ilower,iupper,ictr,local_size
        INTEGER(kind=it4b),DIMENSION(:),allocatable:: rows
        REAL(KIND=dp),DIMENSION(:),allocatable:: rhs,xval
        REAL(KIND=dp):: BetaP,BetaM,BetaW,Lamda
        BetaP = 1.d0/(row/Roref)
        BetaM = 1.d0/(roa/Roref)
        ilower = 0
        iupper = TCell%ExtCell
        local_size = iupper-ilower+1 ! the number of rows
        ! In here, we apply boundary condition for deltaP with its values is 0 at
        ! all boundary. therefore, we do not need to set boundary in vector b
        allocate(rhs(0:TCell%ExtCell))
        allocate(xval(0:TCell%ExtCell))
        allocate(rows(0:TCell%ExtCell))
        call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,b,ierr)
        call HYPRE_IJVectorSetObjectTYPE(b,HYPRE_PARCSR,ierr)
        call HYPRE_IJVectorInitialize(b,ierr)
        call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,x,ierr)
        call HYPRE_IJVectorSetObjectTYPE(x,HYPRE_PARCSR,ierr)
        call HYPRE_IJVectorInitialize(x,ierr)
        rhs(:) = 0.d0
        do i = ibeg,Isize-iu
          do j = jbeg,Jsize-iv
            ictr = TCell%PosNu(i,j)
            rhm(i,j)=0.d0
            if(ictr/=-1) then
              if(iu==1) then
               ! rhs(ictr)=u(i,j)*TGrid%dx(i,j)*TGrid%dy(i,j)/dt     ! Contribution of velocity at old time-step
                rhs(ictr)=rhs(ictr)-IJFlux(i,j)                     ! Convective flux
              else
               ! rhs(ictr)=v(i,j)*TGrid%dx(i,j)*TGrid%dy(i,j)/dt     ! Contribution of velocity at old time-step
                rhs(ictr)=rhs(ictr)-IJFlux(i,j)                     ! Contribution of convective term
           !   Contribution of gravity (should test again 05/10/2017 Contribution of volume flux)
                rhs(ictr)=rhs(ictr)-1.d0/Fr**2.d0*TGrid%dx(i,j)*TGrid%dy(i,j)
              end if
              if(i==ibeg)rhs(ictr)=rhs(ictr)+CWE(j,1)*(dble(iu)*u(i-1,j)+        &
                                                     dble(iv)*v(i-1,j))
              if(i==Isize-iu)rhs(ictr)=rhs(ictr)+CWE(j,2)*(dble(iu)*u(i+1,j)+    &
                                                         dble(iv)*v(i+1,j))
              if(j==jbeg)rhs(ictr)=rhs(ictr)+CSN(i,1)*(dble(iu)*u(i,j-1)+        &
                                                     dble(iv)*v(i,j-1))
              if(j==Jsize-iv)rhs(ictr)=rhs(ictr)+CSN(j,2)*(dble(iu)*u(i,j+1)+    &
                                                         dble(iv)*v(i,j+1))
              if(isnan(rhs(ictr)).or.dabs(rhs(ictr))>1.d20) then
                print*,i,j
                print*,rhs(ictr)
                pause 'Predictor 426'
              end if
              rhm(i,j)=rhs(ictr)
              xval(ictr) = 0.d0
              rows(ictr) = ilower+ictr
            end if
          end do
        end do
        call HYPRE_IJVectorSetValues(b,local_size,rows,rhs,ierr)
        call HYPRE_IJVectorSetValues(x,local_size,rows,xval,ierr)
        call HYPRE_IJVectorAssemble(b,ierr)
        call HYPRE_IJVectorAssemble(x,ierr)
      ! get the x and b objects
        call HYPRE_IJVectorGetObject(b,par_b,ierr)
        call HYPRE_IJVectorGetObject(x,par_x,ierr)
        deallocate(rhs)
        deallocate(xval)
        deallocate(rows)
    end subroutine SetVectors

    SUBROUTINE DeltaGetValues(x,TCell,Var,iu,iv)
        INTEGER*8,INTENT(IN):: x
        INTEGER,INTENT(IN):: iu,iv
        TYPE(Cell),INTENT(IN):: TCell
        REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(INOUT):: Var
        INTEGER(kind=it4b):: i,j
        INTEGER(kind=it4b):: ilower,iupper,local_size,ctr
        INTEGER(kind=it4b),DIMENSION(:),allocatable:: rows
        REAL(KIND=dp),DIMENSION(:),allocatable:: values
        ilower = 0
        iupper = TCell%ExtCell
        local_size = TCell%ExtCell+1 ! number of element
        allocate(values(ilower:iupper))
        allocate(rows(ilower:iupper))
        do i = 1,Isize-iu
          do j = 1,Jsize-iv
            if(TCell%Posnu(i,j)/=-1) then
              rows(TCell%PosNu(i,j)) = TCell%PosNu(i,j)+ilower
            end if
          end do
        end do
        call HYPRE_IJVectorGetValues(x,local_size,rows,values,ierr)
        ctr = 0
        do i = 1,Isize
          do j = 1,Jsize
            if(TCell%PosNu(i,j)==ctr) then
              Var(i,j) = values(ctr)
              ctr = ctr+1
            else
              Var(i,j) = 0.d0
            end if
            if(isnan(Var(i,j))) then
              print*,i,j
              pause 'PredictorUV_564'
            end if
          end do
        end do
        deallocate(values,rows)
    END SUBROUTINE DeltaGetValues

  ! Face flux with MUSCL scheme
    SUBROUTINE ModifiedConvectiveFlux(PGrid,UGrid,VGrid,PCell,UCell,VCell,     &
                                      flux,vb,idir,jdir)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      REAL(KIND=dp),INTENT(IN):: vb
      INTEGER(kind=it4b),INTENT(IN):: idir,jdir
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT):: flux
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp):: uw,vw,us,vs,epsi,uwn,uwp,vsn,vsp
      REAL(KIND=dp):: eta,Sx,Sy
      epsi = 1.d-20
      do i = 1,Isize+idir
        do j = 1,Jsize+jdir
          if(idir==1) then
            uw = 0.5d0*(u(i,j)+u(i-1,j))
     !      For UCell both for convective velocity and scalar velocity
            if(i==ibeg) then
              Flux(i,j,1)=uw**2.d0*UGrid%dy(i,j)*UCell%WEdge_Area(i,j)
            elseif(i>=ibeg+Isize-1) then
              Flux(i,j,1)=uw**2.d0*UGrid%dy(i-1,j)*UCell%EEdge_Area(i-1,j)
            else
              eta=UCell%EtaE(i-1,j)
              uw=(1.d0-eta)*u(i-1,j)+eta*u(i,j)
           !  Second order
           !   Flux(i,j,1)=uw**2.d0*UCell%AlE(i-1,j)**2.d0*                     &
           !                                 UCell%WEdge_Area(i,j)*UGrid%dy(i,j)
           !  First order upwind
              uwp=0.5d0*(uw+dabs(uw))
              uwn=0.5d0*(uw-dabs(uw))
              Flux(i,j,1)=(uwp*u(i-1,j)+uwn*u(i,j))*UCell%AlE(i-1,j)*          &
                                            UGrid%dy(i,j)*UCell%WEdge_Area(i,j)
            end if
            ! Convective velocity: u, scalar advective : v
            uw = 0.5d0*(u(i-1,j+1)+u(i-1,j))
            if(i>ibeg+Isize-1) then
              Flux(i,j,2)=uw*0.5d0*(v(i-1,j)+v(i,j))*VGrid%dy(i-1,j)*          &
                                                     VCell%EEdge_Area(i-1,j)
            elseif(i==ibeg) then
              Flux(i,j,2)=uw*0.5d0*(v(i-1,j)+v(i,j))*VGrid%dy(i,j)*            &
                                                     VCell%WEdge_Area(i,j)
            else
              Sy=UCell%SyN(i-1,j)
              eta=dabs(VCell%FCE(i-1,j,2)+PGrid%dy(i-1,j)/2.d0-                &
                                                   UCell%Cell_Cent(i-1,j,2))/Sy
              if(dabs(eta)>=1.d0) eta=0.5d0
              uw=(1.d0-eta)*u(i-1,j)+eta*u(i-1,j+1)
           !  second order
           !   vw=(1.d0-VCell%EtaE(i-1,j))*v(i-1,j)+VCell%EtaE(i-1,j)*v(i,j)
           !   Flux(i,j,2)=uw*vw*VCell%AlE(i-1,j)*VCell%WEdge_Area(i,j)*        &
           !                                            VGrid%dy(i,j)
           !  firt order upwind
              uwp=0.5d0*(uw+dabs(uw))
              uwn=0.5d0*(uw-dabs(uw))
              Flux(i,j,2)=(uwp*v(i-1,j)+uwn*v(i,j))*VCell%WEdge_Area(i,j)*     &
                                                        VGrid%dy(i,j)
            end if
          end if
          if(jdir==1) then ! Jflux
        ! Convective velocity: v, scalar advective: u
            vs=0.5d0*(v(i,j-1)+v(i+1,j-1))
            if(j>jbeg+Jsize-1) then
              Flux(i,j,1)=vs*0.5d0*(u(i,j)+u(i,j-1))*UGrid%dx(i,j-1)*          &
                                                     UCell%NEdge_Area(i,j-1)
            elseif(j==jbeg) then
              Flux(i,j,1)=vs*0.5d0*(u(i,j)+u(i,j-1))*UGrid%dx(i,j)*            &
                                                     UCell%SEdge_Area(i,j)
            else
              Sx=VCell%SxE(i,j-1)
              eta=dabs(UCell%FCN(i,j-1,1)+PGrid%dx(i,j-1)/2.d0-                &
                                                   VCell%Cell_Cent(i,j-1,1))/Sx
              if(dabs(eta)>=1.d0) eta=0.5d0
              vs=(1.d0-eta)*v(i,j-1)+eta*v(i+1,j-1)
           !  second order
           !   us=(1.d0-UCell%EtaN(i,j-1))*u(i,j-1)+UCell%EtaN(i,j-1)*u(i,j)
           !   Flux(i,j,2)=us*vs*UCell%AlN(i,j-1)*UCell%SEdge_Area(i,j)*        &
           !                                                       UGrid%dx(i,j)
           !  first order upwind
              vsp=0.5d0*(vs+dabs(vs))
              vsn=0.5d0*(vs-dabs(vs))
              Flux(i,j,1)=(vsp*u(i,j-1)+vsn*u(i,j))*UCell%SEdge_Area(i,j)*     &
                                                                 UGrid%dx(i,j)
            end if
            vs = 0.5d0*(v(i,j)+v(i,j-1))
            if(j>=jbeg+Jsize-1) then
              Flux(i,j,2)=vs**2.d0*VGrid%dx(i,Jsize)*VCell%NEdge_Area(i,j-1)
            elseif(j==jbeg) then
              Flux(i,j,2)=vs**2.d0*VGrid%dx(i,j)*VCell%SEdge_Area(i,j)
            else
              eta=VCell%EtaN(i,j-1)
              vs=(1.d0-eta)*v(i,j-1)+eta*v(i,j)
              vs=((vs-vb)*VCell%AlN(i,j-1)+vb)
           !  Second order
           !   Flux(i,j,2)=vs**2.d0*VCell%AlN(i,j-1)**2.d0*VCell%SEdge_Area(i,j)&
           !                 *VGrid%dx(i,j)
           !  First order upwind
              vsp=0.5d0*(vs+dabs(vs))
              vsn=0.5d0*(vs-dabs(vs))
              Flux(i,j,2)=(vsp*v(i,j-1)+vsn*v(i,j))*VCell%SEdge_Area(i,j)*     &
                                                                  VGrid%dx(i,j)
              if(dabs(Flux(i,j,2))>1.d10) then
                print*,vs,eta,vb
                print*,VCell%AlN(i,j-1),VCell%EtaN(i,j-1)
                print*,
              end if
            end if
          end if
        end do
      end do
    END SUBROUTINE ModifiedConvectiveFlux

    SUBROUTINE SecondOrderConvectiveFlux(PGrid,UGrid,VGrid,PCell,UCell,VCell,  &
                                      flux,vb,idir,jdir)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      REAL(KIND=dp),INTENT(IN):: vb
      INTEGER(kind=it4b),INTENT(IN):: idir,jdir
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT):: flux
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp):: xew,xe,xw,yns,yn,ys,vsp,vsn,uwp,uwn
      REAL(KIND=dp):: uw,us,vw,vs,epsil,delh,delh1,delh2,delhec,eta,nx,ny,sx,sy
      epsil = 1.d-20
      do i=1,Isize+idir
        do j=1,Jsize+jdir
          if(idir==1) then
            uw=0.5d0*(u(i,j)+u(i-1,j))
            ! For UCell both for convective velocity and scalar velocity
            if(i==1) then
              Flux(i,j,1)=uw**2.d0*UCell%WEdge_Area(i,j)*UGrid%dy(i,j)
            elseif(i>Isize-1) then
              Flux(i,j,1)=uw**2.d0*UCell%EEdge_Area(i-1,j)*UGrid%dy(i-1,j)
            else
              Sx=UCell%SxE(i-1,j)
              Sy=UCell%Cell_Cent(i,j,2)-UCell%Cell_Cent(i-1,j,2)
              eta=UCell%EtaE(i-1,j)
              uw=(1-eta)*u(i-1,j)+eta*u(i,j)
              Flux(i,j,1)=(uw*UCell%AlE(i-1,j))**2.d0*UCell%WEdge_Area(i,j)*   &
                                                                UGrid%dy(i,j)
              if(i>2.and.i<Isize-1) then
                if(UCell%vofS(i-2,j)>epsi.or.UCell%vofs(i-1,j)>epsi.or.        &
                   UCell%VofS(i,j)>epsi.or.UCell%vofS(i+1,j)>epsi.or.          &
                   UCell%VofS(i+2,j)>epsi) then
            !  first order for x-direction
            !  if(UCell%vofS(i-1,j)>epsi.or.UCell%VofS(i,j)>epsi) then
                  uwp=0.5d0*(uw+dabs(uw))
                  uwn=0.5d0*(uw-dabs(uw))
                  Flux(i,j,1)=(uwp*u(i-1,j)+uwn*u(i,j))*UCell%AlE(i-1,j)*      &
                                         UGrid%dy(i,j)*UCell%WEdge_Area(i,j)
                end if
              end if
            end if
            !  Convective velocity: u, scalar advective : v
            uw=0.5d0*(u(i-1,j+1)+u(i-1,j))
            vw=0.5d0*(v(i-1,j)+v(i,j))
            if(i>Isize) then
              Flux(i,j,2)=uw*vw*VCell%EEdge_Area(i-1,j)*VGrid%dy(i-1,j)
            elseif(i==ibeg) then
              Flux(i,j,2)=uw*vw*VCell%WEdge_Area(i,j)*VGrid%dy(i,j)
            else
              vw=(1.d0-VCell%EtaE(i-1,j))*v(i-1,j)+VCell%EtaE(i-1,j)*v(i,j)
              if(VCell%WEdge_Area(i,j)<0.5d0) then
                delhec=dabs(VCell%FCE(i-1,j,1)*VCell%nxS(i-1,j)+               &
                        VCell%FCE(i-1,j,2)*VCell%nyS(i-1,j)+VCell%phiS(i-1,j))
                if(UCell%MoExCell(i-1,j+1)/=1.and.UCell%VofS(i-1,j+1)<1.d0-epsi)then
                  delh=dabs(UCell%Cell_Cent(i-1,j+1,1)*UCell%nxS(i-1,j+1)+     &
                        UCell%Cell_Cent(i-1,j+1,2)*UCell%nyS(i-1,j+1)+         &
                        UCell%phiS(i-1,j+1))
                  uw=u(i-1,j+1)
                elseif(UCell%MoExCell(i-1,j)/=1.and.UCell%VofS(i-1,j)<1.d0-epsi)then
                  delh=dabs(UCell%Cell_Cent(i-1,j,1)*UCell%nxS(i-1,j)+         &
                        UCell%Cell_Cent(i-1,j,2)*UCell%nyS(i-1,j)+             &
                        UCell%phiS(i-1,j))
                  uw=u(i-1,j)
                else
                  delh=delhec
                end if
             !  if(delhec/delh>1.d0) delh=1.d0/1.d0*delhec
                uw=uw*delhec/delh
             !  Flux(i,j,2)=uw*vw*VCell%AlE(i-1,j)*                            &
             !                                VCell%WEdge_Area(i,j)*VGrid%dy(i,j)
             !  First order of accuracies
                uwp=0.5d0*(uw+dabs(uw))
                uwn=0.5d0*(uw-dabs(uw))
                Flux(i,j,2)=(uwp*v(i-1,j)+uwn*v(i,j))*VCell%WEdge_Area(i,j)*   &
                                                          VGrid%dy(i,j)
              else
                Sy=UCell%SyN(i-1,j)
                eta=dabs(VCell%FCE(i-1,j,2)+VGrid%dy(i-1,j)/2.d0-              &
                                               UCell%Cell_Cent(i-1,j,2))/Sy
                if(dabs(eta)>=1.d0) eta=0.5d0
                uw=(1.d0-eta)*u(i-1,j)+eta*u(i-1,j+1)
                Flux(i,j,2)=uw*vw*VCell%AlE(i-1,j)*VCell%WEdge_Area(i,j)*      &
                                                                 VGrid%dy(i,j)
              ! First order of accuracies
                if(i>2.and.i<Isize-1) then
                  if(VCell%vofS(i-2,j)>epsi.or.VCell%vofs(i-1,j)>epsi.or.      &
                     VCell%VofS(i,j)>epsi.or.VCell%vofS(i+1,j)>epsi.or.        &
                     VCell%VofS(i+2,j)>epsi) then
                    uwp=0.5d0*(uw+dabs(uw))
                    uwn=0.5d0*(uw-dabs(uw))
                    Flux(i,j,2)=(uwp*v(i-1,j)+uwn*v(i,j))*                     &
                                 VCell%WEdge_Area(i,j)* VGrid%dy(i,j)
                  end if
                end if
              end if
            end if
          end if
          if(jdir==1) then ! Jflux
            ! Convective velocity: v, scalar advective: u
            vs=0.5d0*(v(i,j-1)+v(i+1,j-1))
            us=0.5d0*(u(i,j-1)+u(i,j))
            if(j>Jsize) then
              Flux(i,j,1)=vs*us*UCell%NEdge_Area(i,j-1)*UGrid%dx(i,j-1)
            elseif(j==1) then
              Flux(i,j,1)=vs*us*UCell%SEdge_Area(i,j)*UGrid%dx(i,j)
            else
              us=(1.d0-UCell%EtaN(i,j-1))*u(i,j-1)+UCell%EtaN(i,j-1)*u(i,j)
              if(UCell%SEdge_Area(i,j)<0.5d0) then
                delhec=dabs(UCell%FCN(i,j-1,1)*UCell%nxS(i,j-1)+               &
                    UCell%FCN(i,j-1,2)*UCell%nyS(i,j-1)+UCell%phiS(i,j-1))
                if(VCell%MoExCell(i,j-1)/=1.and. VCell%VofS(i,j-1)<1.d0-epsi)then
                  delh=dabs(VCell%Cell_Cent(i,j-1,1)*VCell%nxS(i,j-1)+         &
                    VCell%Cell_Cent(i,j-1,2)*VCell%nyS(i,j-1)+VCell%phiS(i,j-1))
                  vs=v(i,j-1)
                elseif(VCell%MoExCell(i+1,j-1)/=1.and.                         &
                                              VCell%VofS(i+1,j-1)<1.d0-epsi) then
                  delh=dabs(VCell%Cell_Cent(i+1,j-1,1)*VCell%nxS(i+1,j-1)+     &
                    VCell%Cell_Cent(i+1,j-1,2)*VCell%nyS(i+1,j-1)+             &
                    VCell%phiS(i+1,j-1))
                  vs=v(i+1,j-1)
                else
                  delh=delhec
                  vs=vb
                end if
     !           if(delhec/delh>1.d0) delh=1.d0/1.d0*delhec
                vs=vb+(vs-vb)*delhec/(delh+tol)
              !  Flux(i,j,1)=vs*us*UCell%AlN(i,j-1)*                &
              !                              UCell%SEdge_Area(i,j)*UGrid%dx(i,j)
              ! first order of accuracy
                vsp=0.5d0*(vs+dabs(vs))
                vsn=0.5d0*(vs-dabs(vs))
                Flux(i,j,1)=(vsp*u(i,j-1)+vsn*u(i,j))*UCell%SEdge_Area(i,j)*   &
                                                                 UGrid%dx(i,j)
              else
                Sx=VCell%SxE(i,j-1)
                eta=dabs(UCell%FCN(i,j-1,1)+0.5d0*UGrid%dy(i,j-1)-             &
                                              VCell%Cell_Cent(i,j-1,1))/Sx
                if(dabs(eta)>=1.d0) eta=0.5d0
                vs=(1.d0-eta)*v(i,j-1)+eta*v(i+1,j-1)
                Flux(i,j,1)=vs*us*UCell%AlN(i,j-1)*UCell%SEdge_Area(i,j)*      &
                                                                 UGrid%dx(i,j)
                if(j>2.and.j<Jsize-1) then
                  if(UCell%vofS(i,j-2)>epsi.or.UCell%vofs(i,j-1)>epsi.or.      &
                     UCell%VofS(i,j)>epsi.or.UCell%vofS(i,j+1)>epsi.or.        &
                     UCell%vofS(i,j+2)>epsi) then
              ! First order in y direction
                    vsp=0.5d0*(vs+dabs(vs))
                    vsn=0.5d0*(vs-dabs(vs))
                    Flux(i,j,1)=(vsp*u(i,j-1)+vsn*u(i,j))*                     &
                                 UCell%SEdge_Area(i,j)*UGrid%dx(i,j)
                  end if
                endif
              end if
            end if
            vs=0.5d0*(v(i,j)+v(i,j-1))
            if(j>Jsize-1) then
              Flux(i,j,2)=vs**2.d0*VCell%NEdge_Area(i,j-1)*VGrid%dx(i,j-1)
            elseif(j==1) then
              Flux(i,j,2)=vs**2.d0*VCell%SEdge_Area(i,j)*VGrid%dx(i,j)
            else
              Sx=VCell%Cell_Cent(i,j,1)-VCell%Cell_Cent(i,j-1,1)
              Sy=VCell%SyN(i,j-1)
              eta=VCell%EtaN(i,j-1)
              vs=(1.d0-eta)*v(i,j-1)+eta*v(i,j)
              vs=((vs-vb)*VCell%AlN(i,j-1)+vb)
              Flux(i,j,2)=vs**2.d0*VCell%SEdge_Area(i,j)*VGrid%dx(i,j)
              if(j>2.and.j<Jsize-1) then
                if(VCell%vofS(i,j-2)>epsi.or.VCell%vofs(i,j-1)>epsi.or.        &
                   VCell%VofS(i,j)>epsi.or.VCell%vofS(i,j+1)>epsi.or.          &
                   VCell%vofS(i,j+2)>epsi) then
          !   First order of accuracy
                  vsp=0.5d0*(vs+dabs(vs))
                  vsn=0.5d0*(vs-dabs(vs))
                  Flux(i,j,2)=(vsp*v(i,j-1)+vsn*v(i,j))*VCell%SEdge_Area(i,j)* &
                                                      VGrid%dx(i,j)
                end if
              end if
            end if
          end if
        end do
      end do
    END SUBROUTINE SecondOrderConvectiveFlux

    SUBROUTINE HighOrderConvectiveFluxForXDir(PGrid,UGrid,VGrid,UCell,         &
                                      VCell,vb,flux)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: UCell,VCell
      REAL(KIND=dp),INTENT(IN):: vb
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT)::flux
      INTEGER(kind=it4b):: i,j,Lim
      REAL(KIND=dp):: ul,ur,vl,vr,alr,uwp,uwn,sx,sy,uw,delhec,delh,eta
      REAL(KIND=dp):: omei,omei1,ri,ri1,tolim,tol
      tol=1.d-24
      Lim=2
      do j=1,Jsize
        ul=u(0,j);ur=u(1,j)
        alr=(ur+ul)
        alr=dmax1(dabs(ur),dabs(ul))
        flux(1,j,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*               &
                                             UCell%WEdge_Area(1,j)*UGrid%dy(1,j)
        vl=v(0,j);vr=v(1,j)
        alr=0.5d0*(u(0,j)+u(0,j+1))
        flux(1,j,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*                   &
                                             VCell%WEdge_Area(1,j)*VGrid%dy(1,j)
        do i=2,Isize
        ! Calculate threshold for MUSCL
        ! from 'A MUSCL scheme on staggered grids with kinetic-like fluxes
        ! for the barotropic Euler system', Thierry Goundon, Julie Llobell
          tolim=dmin1(2.d0*UGrid%dx(i-1,j)/PGrid%dx(i,j),                      &
                      2.d0*UGrid%dx(i,j)/PGrid%dx(i,j))
          ri=((u(i-1,j)-u(i-2,j))/PGrid%dx(i-1,j))/                            &
             ((u(i,j)-u(i-1,j))/PGrid%dx(i,j))
          omei=MUSCLLimiter(ri,Lim,tolim)*(u(i,j)-u(i-1,j))/PGrid%dx(i,j)
          ri1=((u(i,j)-u(i-1,j))/PGrid%dx(i,j))/                               &
              ((u(i+1,j)-u(i,j))/PGrid%dx(min(Isize,i+1),j))
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(u(i+1,j)-u(i,j))/                 &
                                                      PGrid%dx(min(Isize,i+1),j)
          ul=u(i-1,j)+0.5d0*(PGrid%dx(i,j))*omei
          ur=u(i,j)-0.5d0*(PGrid%dx(i,j))*omei1
          alr=(ur+ul)
          alr=dmax1(dabs(ur),dabs(ul))
          flux(i,j,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*             &
                                             UCell%WEdge_Area(i,j)*UGrid%dy(i,j)
          if(i>2.and.i<Isize) then
            if(UCell%vofS(i-2,j)>epsi.or.UCell%vofs(i-1,j)>epsi.or.            &
                           UCell%VofS(i,j)>epsi.or.UCell%vofS(i+1,j)>epsi) then
              Sx=UCell%SxE(i-1,j)
              Sy=UCell%Cell_Cent(i,j,2)-UCell%Cell_Cent(i-1,j,2)
              eta=UCell%EtaE(i-1,j)
              uw=(1-eta)*u(i-1,j)+eta*u(i,j)
        !   first order for x-direction
        !   if(UCell%vofS(i-1,j)>epsi.or.UCell%VofS(i,j)>epsi) then
              uwp=0.5d0*(uw+dabs(uw))
              uwn=0.5d0*(uw-dabs(uw))
              Flux(i,j,1)=(uwp*u(i-1,j)+uwn*u(i,j))*UCell%AlE(i-1,j)*            &
                                        UGrid%dy(i,j)*UCell%WEdge_Area(i,j)
            end if
          end if
          tolim=dmin1(2.d0*VGrid%dx(i-1,j)/UGrid%dx(i-1,j),                    &
                      2.d0*VGrid%dx(i,j)/UGrid%dx(i-1,j))
          ri=((v(i-1,j)-v(i-2,j))/UGrid%dx(max(1,i-2),j))/                     &
             ((v(i,j)-v(i-1,j))/UGrid%dx(i-1,j))
          omei=MUSCLlimiter(ri,Lim,tolim)*(v(i,j)-v(i-1,j))/UGrid%dx(i-1,j)
          ri1=((v(i,j)-v(i-1,j))/UGrid%dx(i-1,j))/                             &
              ((v(i+1,j)-v(i,j))/UGrid%dx(i,j))
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(v(i+1,j)-v(i,j))/UGrid%dx(i,j)
          vl=v(i-1,j)+0.5d0*(PGrid%dx(i-1,j))*omei
          vr=v(i,j)-0.50d0*(PGrid%dx(i,j))*omei1
          alr=0.5d0*(u(i-1,j)+u(i-1,j+1))
        !  alr=dmax1(dabs(u(i-1,j)),dabs(u(i-1,j+1)))
          flux(i,j,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*                 &
                                             VCell%WEdge_Area(i,j)*VGrid%dy(i,j)
          if(i>2.and.i<Isize) then
            if(VCell%vofS(i-2,j)>epsi.or.VCell%vofs(i-1,j)>epsi.or.            &
                           VCell%VofS(i,j)>epsi.or.VCell%vofS(i+1,j)>epsi) then
              if(VCell%WEdge_Area(i,j)<0.5d0) then
                delhec=dabs(VCell%FCE(i-1,j,1)*VCell%nxS(i-1,j)+                 &
                      VCell%FCE(i-1,j,2)*VCell%nyS(i-1,j)+VCell%phiS(i-1,j))
                if(UCell%MoExCell(i-1,j+1)/=1.and.UCell%VofS(i-1,j+1)<1.d0-epsi)then
                  delh=dabs(UCell%Cell_Cent(i-1,j+1,1)*UCell%nxS(i-1,j+1)+       &
                      UCell%Cell_Cent(i-1,j+1,2)*UCell%nyS(i-1,j+1)+           &
                      UCell%phiS(i-1,j+1))+tol
                  uw=u(i-1,j+1)
                elseif(UCell%MoExCell(i-1,j)/=1.and.UCell%VofS(i-1,j)<1.d0-epsi)then
                  delh=dabs(UCell%Cell_Cent(i-1,j,1)*UCell%nxS(i-1,j)+           &
                      UCell%Cell_Cent(i-1,j,2)*UCell%nyS(i-1,j)+               &
                      UCell%phiS(i-1,j))+tol
                  uw=u(i-1,j)
                else
                  delh=delhec+tol
                  uw=0.d0
                end if
                uwp=0.5d0*(uw+dabs(uw))
                uwn=0.5d0*(uw-dabs(uw))
                Flux(i,j,2)=(uwp*v(i-1,j)+uwn*v(i,j))*delhec/delh*               &
                                             VCell%WEdge_Area(i,j)*VGrid%dy(i,j)
              else
                Sy=UCell%SyN(i-1,j)
                eta=dabs(VCell%FCE(i-1,j,2)+VGrid%dy(i-1,j)/2.d0-                &
                                               UCell%Cell_Cent(i-1,j,2))/Sy
                if(dabs(eta)>=1.d0) eta=0.5d0
                uw=(1.d0-eta)*u(i-1,j)+eta*u(i-1,j+1)
                uwp=0.5d0*(uw+dabs(uw))
                uwn=0.5d0*(uw-dabs(uw))
                Flux(i,j,2)=(uwp*v(i-1,j)+uwn*v(i,j))*VCell%WEdge_Area(i,j)*     &
                                                          VGrid%dy(i,j)

              end if
            end if
          end if
        end do
        ul=u(Isize,j);ur=u(Isize+1,j)
        alr=(ur+ul)
        alr=dmax1(dabs(ur),dabs(ul))
        flux(Isize+1,j,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*         &
                                     UCell%EEdge_Area(Isize,j)*UGrid%dy(Isize,j)
        vl=v(Isize,j);vr=v(Isize+1,j)
        alr=0.5d0*(u(Isize,j)+u(Isize,j+1))
        flux(Isize+1,j,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*             &
                                     VCell%EEdge_Area(Isize,j)*VGrid%dy(Isize,j)
      end do
    end subroutine HighOrderConvectiveFluxForXDir

    subroutine HighOrderConvectiveFluxForYDir(PGrid,UGrid,VGrid,UCell,         &
                                      VCell,vb,flux)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: UCell,VCell
      REAL(KIND=dp),INTENT(IN):: vb
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT)::flux
      INTEGER(kind=it4b):: i,j,Lim
      REAL(KIND=dp):: ul,ur,vl,vr,alr,delhec,delh,tol
      REAL(KIND=dp):: omei,omei1,ri,ri1,tolim,Sx,Sy,vs,vsp,vsn,eta
      Lim=2
      tol=1.d-24
      do i=1,Isize
        ul=u(i,0);ur=u(i,1)
        alr=0.5d0*(v(i,0)+v(i+1,0))
        flux(i,1,1)=0.5d0*(alr*ur*+alr*ul-dabs(alr)*(ur-ul))*                  &
                                            UCell%SEdge_Area(i,1)*UGrid%dx(i,1)
        vl=v(i,0);vr=v(i,1)
        alr=(vr+vl)
        alr=dmax1(dabs(vl),dabs(vr))
        flux(i,1,2)=0.5d0*(vr**2.d0+vl**2.d0-dabs(alr)*(vr-vl))*               &
                                            VCell%SEdge_Area(i,1)*VGrid%dx(i,1)
        do j=2,Jsize
          tolim=dmin1(2.d0*UGrid%dy(i,j-1)/VGrid%dy(i,j-1),                    &
                      2.d0*UGrid%dy(i,j)/VGrid%dy(i,j-1))
          ri=((u(i,j-1)-u(i,j-2))/VGrid%dy(i,max(1,j-2)))/                     &
             ((u(i,j)-u(i,j-1))/VGrid%dy(i,j-1))
          omei=MUSCLLimiter(ri,Lim,tolim)*(u(i,j)-u(i,j-1))/VGrid%dy(i,j-1)
          ri1=((u(i,j)-u(i,j-1))/VGrid%dy(i,j-1))/                             &
              ((u(i,j+1)-u(i,j))/VGrid%dy(i,j))
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(u(i,j+1)-u(i,j))/VGrid%dy(i,j)
          ul=u(i,j-1)+0.5d0*(PGrid%dy(i,j-1))*omei
          ur=u(i,j)-0.5d0*(PGrid%dy(i,j))*omei1
          alr=0.5d0*(v(i,j-1)+v(i+1,j-1))
       !   alr=dmax1(dabs(v(i,j-1)),dabs(v(i+1,j-1)))
          flux(i,j,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*                 &
                                             UCell%SEdge_Area(i,j)*UGrid%dx(i,j)
          if(j>2.and.j<Jsize) then
            if(UCell%vofS(i,j-2)>epsi.or.UCell%vofs(i,j-1)>epsi.or.            &
                           UCell%VofS(i,j)>epsi.or.UCell%vofS(i,j+1)>epsi) then
              if(UCell%SEdge_Area(i,j)<0.5d0) then
                delhec=dabs(UCell%FCN(i,j-1,1)*UCell%nxS(i,j-1)+               &
                   UCell%FCN(i,j-1,2)*UCell%nyS(i,j-1)+UCell%phiS(i,j-1))
                if(VCell%MoExCell(i,j-1)/=1.and. VCell%VofS(i,j-1)<1.d0-epsi)then
                  delh=dabs(VCell%Cell_Cent(i,j-1,1)*VCell%nxS(i,j-1)+         &
                    VCell%Cell_Cent(i,j-1,2)*VCell%nyS(i,j-1)+                 &
                    VCell%phiS(i,j-1))+tol
                  vs=v(i,j-1)
                elseif(VCell%MoExCell(i+1,j-1)/=1.and.                         &
                                            VCell%VofS(i+1,j-1)<1.d0-epsi) then
                  delh=dabs(VCell%Cell_Cent(i+1,j-1,1)*VCell%nxS(i+1,j-1)+     &
                  VCell%Cell_Cent(i+1,j-1,2)*VCell%nyS(i+1,j-1)+               &
                  VCell%phiS(i+1,j-1))+tol
                  vs=v(i+1,j-1)
                else
                  delh=delhec+tol
                  vs=0.d0
                end if
                vs=vb+(vs-vb)*delhec/delh
                vsp=0.5d0*(vs+dabs(vs))
                vsn=0.5d0*(vs-dabs(vs))
                Flux(i,j,1)=(vsp*u(i,j-1)+vsn*u(i,j))*delhec/delh*               &
                                            UCell%SEdge_Area(i,j)*UGrid%dx(i,j)
                if(isnan(flux(i,j,1))) then
                  print*,vsp,vsn
                  print*,delhec,delh,UCell%SEdge_Area(i,j)
                  print*,'fuck you bugs'
                end if
              else
                Sx=VCell%SxE(i,j-1)
                eta=dabs(UCell%FCN(i,j-1,1)+0.5d0*UGrid%dy(i,j-1)-               &
                                            VCell%Cell_Cent(i,j-1,1))/Sx
                if(dabs(eta)>=1.d0) eta=0.5d0
                vs=(1.d0-eta)*v(i,j-1)+eta*v(i+1,j-1)
                vsp=0.5d0*(vs+dabs(vs))
                vsn=0.5d0*(vs-dabs(vs))
                Flux(i,j,1)=(vsp*u(i,j-1)+vsn*u(i,j))*UCell%SEdge_Area(i,j)*     &
                                                                 UGrid%dx(i,j)
              end if
            end if
          end if

          tolim=dmin1(2.d0*VGrid%dy(i,j-1)/PGrid%dy(i,j),                      &
                      2.d0*VGrid%dy(i,j)/PGrid%dy(i,j))
          ri=((v(i,j-1)-v(i,j-2))/PGrid%dy(i,j-1))/                            &
             ((v(i,j)-v(i,j-1))/PGrid%dy(i,j))
          omei=MUSCLlimiter(ri,Lim,tolim)*(v(i,j)-v(i,j-1))/PGrid%dy(i,j-1)
          ri1=((v(i,j)-v(i,j-1))/PGrid%dy(i,j))/                               &
              ((v(i,j+1)-v(i,j))/PGrid%dy(i,min(Jsize,j+1)))
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(v(i,j+1)-v(i,j))/                 &
                                               PGrid%dy(i,min(Jsize,j+1))
          vl=v(i,j-1)+0.5d0*(PGrid%dy(i,j))*omei
          vr=v(i,j)-0.5d0*(PGrid%dy(i,j))*omei1
          alr=(vr+vl)
          alr=dmax1(dabs(vl),dabs(vr))
          flux(i,j,2)=0.5d0*(vr*vr+vl*vl-dabs(alr)*(vr-vl))*                   &
                                             VCell%SEdge_Area(i,j)*VGrid%dx(i,j)
          if(j>2.and.j<Jsize) then
            if(VCell%vofS(i,j-2)>epsi.or.VCell%vofs(i,j-1)>epsi.or.            &
                           VCell%VofS(i,j)>epsi.or.VCell%vofS(i,j+1)>epsi) then
              Sx=VCell%Cell_Cent(i,j,1)-VCell%Cell_Cent(i,j-1,1)
              Sy=VCell%SyN(i,j-1)
              eta=VCell%EtaN(i,j-1)
              vs=(1.d0-eta)*v(i,j-1)+eta*v(i,j)
              vs=(vb+(vs-vb)*VCell%AlN(i,j-1))
              vsp=0.5d0*(vs+dabs(vs))
              vsn=0.5d0*(vs-dabs(vs))
              Flux(i,j,2)=(vsp*v(i,j-1)+vsn*v(i,j))*                            &
                                             VCell%SEdge_Area(i,j)*VGrid%dx(i,j)
            end if
          end if
        end do
        ul=u(i,Jsize);ur=u(i,Jsize+1)
        alr=0.5d0*(v(i,Jsize)+v(i+1,Jsize))
        flux(i,Jsize+1,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*             &
                                     UCell%NEdge_Area(i,Jsize)*UGrid%dx(i,Jsize)
        vl=v(i,Jsize);vr=v(i,Jsize+1)
        alr=(vr+vl)
        alr=dmax1(dabs(vl),dabs(vr))
        flux(i,Jsize+1,2)=0.5d0*(vr**2.d0+vl**2.d0-dabs(alr)*(vr-vl))*         &
                                     VCell%NEdge_Area(i,Jsize)*VGrid%dx(i,Jsize)
      end do
    end subroutine HighOrderConvectiveFluxForYDir

    subroutine HighOrderDensityBasedConvFluxXDir(PGrid,UGrid,VGrid,PCell,UCell,&
                                      VCell,Uro,Vro,flux,Vflux,dt)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(IN):: Uro,Vro
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT):: Vflux,flux
      REAL(KIND=dp),INTENT(IN):: dt
      INTEGER(kind=it4b):: i,j,Lim
      REAL(KIND=dp):: ul,ur,vl,vr,alr,volf,vols,rolr,rol,ror,uw
      REAL(KIND=dp):: omei,omei1,ri,ri1,tolim,tol
      tol=1.d-24
      Lim=2
      do j=1,Jsize
        if(UCell%WEdge_Area(1,j)>=epsi) then
          rolr=((1.d0-PCell%Vof(1,j)-PCell%VofS(1,j))*roa/Roref+               &
                                                     PCell%Vof(1,j)*row/Roref)
        else
          rolr=0.d0
        end if
        ul=u(0,j);ur=u(1,j)
        alr=(ur+ul)
        if(alr>=0.d0) then
          VFlux(1,j,1)=ul*rolr*UGrid%dy(1,j)
        else
          VFlux(1,j,1)=ur*rolr*UGrid%dy(1,j)
        end if
        flux(1,j,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*               &
                                                             rolr*UGrid%dy(1,j)
        vl=v(0,j);vr=v(1,j)
        alr=0.5d0*(u(0,j)+u(0,j+1))
        if(VCell%WEdge_ARea(1,j)>=epsi) then
          rolr=Vro(1,j)
        else
          rolr=0.d0
        end if
        if(alr>=0.d0) then
          VFlux(1,j,2)=vl*rolr*VGrid%dy(1,j)
        else
          VFlux(1,j,2)=vr*rolr*VGrid%dy(1,j)
        end if
        flux(1,j,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*rolr*VGrid%dy(1,j)
        do i=2,Isize
        ! For UCell
          tolim=dmin1(2.d0*UGrid%dx(i-1,j)/PGrid%dx(i,j),                      &
                      2.d0*UGrid%dx(i,j)/PGrid%dx(i,j))
          ri=((u(i-1,j)-u(i-2,j))/PGrid%dx(i-1,j))/                            &
             ((u(i,j)-u(i-1,j))/PGrid%dx(i,j)+tol)
          omei=MUSCLLimiter(ri,Lim,tolim)*(u(i,j)-u(i-1,j))/PGrid%dx(i,j)
          ri1=((u(i,j)-u(i-1,j))/PGrid%dx(i,j))/                               &
              ((u(i+1,j)-u(i,j))/PGrid%dx(min(Isize,i+1),j)+tol)
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(u(i+1,j)-u(i,j))/                 &
                                                      PGrid%dx(min(Isize,i+1),j)
          ul=u(i-1,j)+0.5d0*(PGrid%dx(i,j))*omei
          ur=u(i,j)-0.5d0*(PGrid%dx(i,j))*omei1
          alr=(ur+ul)
       !  find velocity to drive volume fraction
          if(alr>=0.d0) then
            uw=ul
          else
            uw=ur
          end if
          alr=dmax1(dabs(ur),dabs(ul))
       !  find the density at cell face
          if(UCell%WEdge_Area(i,j)>=epsi) then
            if(PCell%vof(i,j)>1.d0-epsi.or.PCell%vof(i,j)<epsi.or.             &
                                                              dabs(uw)<tol)then
              rolr=((1.d0-PCell%Vof(i,j)-PCell%VofS(i,j))*roa/Roref+           &
                                      PCell%Vof(i,j)*row/Roref)
            else
              if(PCell%VofS(i,j)<epsi) then
                call frac(PCell%nx(i,j),PCell%ny(i,j),PCell%phi(i,j)+          &
                     uw*dt/2.d0*PCell%nx(i,j),dabs(uw*dt),PGrid%dy(i,j),volf)
                rolr=(volf*row/Roref+(1.d0-volf-vols)*roa/Roref)
              else
               call CellGeoCal(PCell%nxs(i,j),PCell%nys(i,j),                 &
                     PCell%phis(i,j)+uw*dt/2.d0*PCell%nxs(i,j),PCell%nx(i,j),  &
                     PCell%ny(i,j),PCell%phi(i,j)+uw*dt/2.d0*PCell%nx(i,j),    &
                                    dabs(uw*dt),PGrid%dy(i,j),vols,volf)
                rolr=(volf*row/Roref+(1.d0-volf-vols)*roa/Roref)
              end if
            end if
      !    no reconstruction need for computing face flux
      !      rolr=((1.d0-PCell%Vof(i,j)-PCell%VofS(i,j))*roa/Roref+           &
      !                                PCell%Vof(i,j)*row/Roref)
          else
            rolr=0.d0
          end if
          if(alr>=0.d0) then
            VFlux(i,j,1)=ul*rolr*UGrid%dy(i,j)
          else
            VFlux(i,j,1)=ur*rolr*UGrid%dy(i,j)
          end if
          Flux(i,j,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*             &
                                                             rolr*UGrid%dy(i,j)
        ! For VCell
          tolim=dmin1(2.d0*VGrid%dx(i-1,j)/UGrid%dx(i-1,j),                    &
                      2.d0*VGrid%dx(i,j)/UGrid%dx(i-1,j))
          ri=((v(i-1,j)-v(i-2,j))/UGrid%dx(max(1,i-2),j))/                     &
             ((v(i,j)-v(i-1,j))/UGrid%dx(i-1,j)+tol)
          omei=MUSCLlimiter(ri,Lim,tolim)*(v(i,j)-v(i-1,j))/UGrid%dx(i-1,j)
          ri1=((v(i,j)-v(i-1,j))/UGrid%dx(i-1,j))/                             &
              ((v(i+1,j)-v(i,j))/UGrid%dx(i,j)+tol)
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(v(i+1,j)-v(i,j))/UGrid%dx(i,j)
          vl=v(i-1,j)+0.5d0*(PGrid%dx(i-1,j))*omei
          vr=v(i,j)-0.50d0*(PGrid%dx(i,j))*omei1
          alr=0.5d0*(u(i-1,j)+u(i-1,j+1))

          ri=((Vro(i-1,j)-Vro(max(1,i-2),j))/UGrid%dx(max(1,i-2),j))/          &
             ((Vro(i,j)-Vro(i-1,j))/UGrid%dx(i-1,j)+tol)
          omei=MUSCLlimiter(ri,Lim,tolim)*(Vro(i,j)-Vro(i-1,j))/UGrid%dx(i-1,j)
          ri1=((Vro(i,j)-Vro(i-1,j))/UGrid%dx(i-1,j))/                         &
              ((Vro(min(Isize,i+1),j)-Vro(i,j))/UGrid%dx(i,j)+tol)
          omei1=MUSCLLimiter(ri1,Lim,tolim)*                                   &
                                 (Vro(min(Isize,i+1),j)-Vro(i,j))/UGrid%dx(i,j)
          rol=Vro(i-1,j)+0.5d0*(PGrid%dx(i-1,j))*omei
          ror=Vro(i,j)-0.50d0*(PGrid%dx(i,j))*omei1
          if(VCell%WEdge_Area(i,j)>=epsi) then
            Vflux(i,j,2)=0.5d0*(alr*ror+alr*rol-dabs(alr)*(ror-rol))*          &
                                                                  VGrid%dy(i,j)
            Flux(i,j,2)=0.50d0*(alr*vr*ror+alr*vl*rol-                         &
                                       dabs(alr)*(ror*vr-rol*vl))*VGrid%dy(i,j)
          else
            Vflux(i,j,2)=0.d0
            Flux(i,j,2)=0.d0
          end if
        end do

      ! For (Isize+1,:) elements of flux
        if(UCell%EEdge_Area(Isize,j)>=epsi) then
          rolr=((1.d0-PCell%Vof(Isize,j)-PCell%VofS(Isize,j))*roa/Roref+       &
                                                 PCell%Vof(Isize,j)*row/Roref)
        else
          rolr=0.d0
        end if
        ul=u(Isize,j);ur=u(Isize+1,j)
        alr=(ur+ul)
        if(alr>=0.d0) then
          VFlux(Isize+1,j,1)=ul*rolr*UGrid%dy(Isize,j)
        else
          VFlux(Isize+1,j,1)=ur*rolr*UGrid%dy(Isize,j)
        end if
        flux(Isize+1,j,1)=0.5d0*(ur**2.d0+ul**2.d0-dabs(alr)*(ur-ul))*         &
                                                         rolr*UGrid%dy(Isize,j)
        vl=v(Isize,j);vr=v(Isize+1,j)
        alr=0.5d0*(u(Isize,j)+u(Isize,j+1))
        if(VCell%EEdge_ARea(Isize,j)>=epsi) then
          rolr=Vro(Isize,j)
        else
          rolr=0.d0
        end if
        if(alr>=0.d0) then
          VFlux(Isize+1,j,2)=vl*rolr*VGrid%dy(Isize,j)
        else
          VFlux(Isize+1,j,2)=vr*rolr*VGrid%dy(Isize,j)
        end if
        flux(Isize+1,j,2)=0.5d0*(alr*vr+alr*vl-dabs(alr)*(vr-vl))*rolr*        &
                                                              VGrid%dy(Isize,j)
      end do
    end subroutine HighOrderDensityBasedConvFluxXDir

    subroutine HighOrderDensityBasedConvFluxYDir(PGrid,UGrid,VGrid,PCell,UCell,&
                                      VCell,Uro,Vro,flux,Vflux,dt)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(IN):: Uro,Vro
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT):: Vflux,flux
      REAL(KIND=dp),INTENT(IN):: dt
      INTEGER(kind=it4b):: i,j,Lim
      REAL(KIND=dp):: ul,ur,vl,vr,alr,rolr,rol,ror,volf,vols,vs
      REAL(KIND=dp):: omei,omei1,ri,ri1,tolim,tol
      tol=1.d-24
      Lim=2
      do i=1,Isize
        if(UCell%SEdge_Area(i,1)>=epsi) then
          rolr=Uro(i,1)
        else
          rolr=0.d0
        end if
        ul=u(i,0);ur=u(i,1)
        alr=0.5d0*(v(i,0)+v(i+1,0))
        if(alr>=0.d0) then
          VFlux(i,1,1)=ul*rolr*UGrid%dx(i,1)
        else
          VFlux(i,1,1)=ur*rolr*UGrid%dx(i,1)
        end if
        flux(i,1,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*rolr*UGrid%dx(i,1)

        vl=v(i,0);vr=v(i,1)
        alr=(vl+vr)
        if(VCell%SEdge_ARea(i,1)>=epsi) then
          rolr=((1.d0-PCell%Vof(i,1)-PCell%VofS(i,1))*roa/Roref+               &
                                                     PCell%Vof(i,1)*row/Roref)
        else
          rolr=0.d0
        end if
        if(alr>=0.d0) then
          VFlux(i,1,2)=vl*rolr*VGrid%dx(i,1)
        else
          VFlux(i,1,2)=vr*rolr*VGrid%dx(i,1)
        end if
        flux(i,1,2)=0.5d0*(vr*vr+vl*vl-dabs(alr)*(vr-vl))*rolr*VGrid%dx(i,1)
        do j=2,Jsize
          tolim=dmin1(2.d0*UGrid%dy(i,j-1)/VGrid%dy(i,j-1),                    &
                      2.d0*UGrid%dy(i,j)/VGrid%dy(i,j-1))
          ri=((u(i,j-1)-u(i,j-2))/VGrid%dy(i,max(1,j-2)))/                     &
             ((u(i,j)-u(i,j-1))/VGrid%dy(i,j-1)+tol)
          omei=MUSCLLimiter(ri,Lim,tolim)*(u(i,j)-u(i,j-1))/VGrid%dy(i,j-1)
          ri1=((u(i,j)-u(i,j-1))/VGrid%dy(i,j-1))/                             &
              ((u(i,j+1)-u(i,j))/VGrid%dy(i,j)+tol)
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(u(i,j+1)-u(i,j))/VGrid%dy(i,j)
          ul=u(i,j-1)+0.5d0*(PGrid%dy(i,j-1))*omei
          ur=u(i,j)-0.5d0*(PGrid%dy(i,j))*omei1
          alr=0.5d0*(v(i,j-1)+v(i+1,j-1))

          ri=((Uro(i,j-1)-Uro(i,max(1,j-2)))/VGrid%dy(i,max(1,j-2)))/          &
             ((Uro(i,j)-Uro(i,j-1))/VGrid%dy(i,j-1)+tol)
          omei=MUSCLLimiter(ri,Lim,tolim)*(Uro(i,j)-Uro(i,j-1))/VGrid%dy(i,j-1)
          ri1=((Uro(i,j)-Uro(i,j-1))/VGrid%dy(i,j-1))/                         &
              ((Uro(i,min(Jsize,j+1))-Uro(i,j))/VGrid%dy(i,j)+tol)
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(Uro(i,min(Jsize,j+1))-Uro(i,j))/  &
                                                               VGrid%dy(i,j)
          rol=Uro(i,j-1)+0.5d0*(PGrid%dy(i,j-1))*omei
          ror=Uro(i,j)-0.5d0*(PGrid%dy(i,j))*omei1
          if(UCell%SEdge_Area(i,j)>=epsi) then
            Vflux(i,j,1)=0.5d0*(ror*alr+rol*alr-dabs(alr)*(ror-rol))*          &
                                                                  UGrid%dx(i,j)
            Flux(i,j,1)=0.5d0*(ror*alr*ur+rol*alr*ul-                          &
                              dabs(alr)*(ror*ur-rol*ul))*UGrid%dx(i,j)
          else
            Vflux(i,j,1)=0.d0
            Flux(i,j,1)=0.d0
          end if
          tolim=dmin1(2.d0*VGrid%dy(i,j-1)/PGrid%dy(i,j),                      &
                      2.d0*VGrid%dy(i,j)/PGrid%dy(i,j))
          ri=((v(i,j-1)-v(i,j-2))/PGrid%dy(i,j-1))/                            &
             ((v(i,j)-v(i,j-1))/PGrid%dy(i,j)+tol)
          omei=MUSCLlimiter(ri,Lim,tolim)*(v(i,j)-v(i,j-1))/PGrid%dy(i,j-1)
          ri1=((v(i,j)-v(i,j-1))/PGrid%dy(i,j))/                               &
              ((v(i,j+1)-v(i,j))/PGrid%dy(i,min(Jsize,j+1))+tol)
          omei1=MUSCLLimiter(ri1,Lim,tolim)*(v(i,j+1)-v(i,j))/                 &
                                               PGrid%dy(i,min(Jsize,j+1))
          vl=v(i,j-1)+0.5d0*(PGrid%dy(i,j))*omei
          vr=v(i,j)-0.5d0*(PGrid%dy(i,j))*omei1
          alr=(vr+vl)

       !  find velocity to drive volume fraction
          if(alr>=0.d0) then
            vs=vl
          else
            vs=vr
          end if
          if(VCell%SEdge_Area(i,j)>=epsi) then
            if(PCell%vof(i,j)>1.d0-epsi.or.PCell%vof(i,j)<epsi.or.             &
                                                              dabs(vs)<tol)then
              rolr=(PCell%vof(i,j)*row/Roref+                                  &
                           (1.d0-PCell%vof(i,j)-PCell%vofS(i,j))*roa/Roref)
            else
              if(PCell%vofS(i,j)<epsi) then
                call frac(PCell%nx(i,j),PCell%ny(i,j),PCell%phi(i,j)+          &
                       vs*dt/2.d0*PCell%ny(i,j),PGrid%dx(i,j),dabs(vs*dt),volf)
                rolr=(volf*row/Roref+(1.d0-volf-vols)*roa/Roref)
              else
                call CellGeoCal(PCell%nxs(i,j),PCell%nys(i,j),                 &
                     PCell%phis(i,j)+vs*dt/2.d0*PCell%nys(i,j),PCell%nx(i,j),  &
                     PCell%ny(i,j),PCell%phi(i,j)+vs*dt/2.d0*PCell%ny(i,j),    &
                                   PGrid%dx(i,j),dabs(vs*dt),vols,volf)
                rolr=(volf*row/Roref+(1.d0-volf-vols)*roa/Roref)
              end if
            end if
       !    no reconstruction need for computing face flux for VCell
       !     rolr=(PCell%vof(i,j)*row/Roref+                                  &
       !                    (1.d0-PCell%vof(i,j)-PCell%vofS(i,j))*roa/Roref)
          else
            rolr=0.d0
          end if
          if(alr>=0.d0) then
            VFlux(i,j,2)=vl*rolr*VGrid%dx(i,j)
          else
            VFlux(i,j,2)=vr*rolr*VGrid%dx(i,j)
          end if
          Flux(i,j,2)=0.5d0*(vr**2.d0+vl**2.d0-dabs(alr)*(vr-vl))*             &
                                                             rolr*VGrid%dx(i,j)
        end do
     !  For (:,Jsize+1) elements of flux
        if(UCell%NEdge_Area(i,Jsize)>=epsi) then
          rolr=Uro(i,Jsize)
        else
          rolr=0.d0
        end if
        ul=u(i,Jsize);ur=u(i,Jsize+1)
        alr=0.5d0*(v(i,Jsize)+v(i+1,Jsize))
        if(alr>=0.d0) then
          VFlux(i,Jsize+1,1)=ul*rolr*UGrid%dx(i,Jsize)
        else
          VFlux(i,Jsize+1,1)=ur*rolr*UGrid%dx(i,Jsize)
        end if
        flux(i,Jsize+1,1)=0.5d0*(alr*ur+alr*ul-dabs(alr)*(ur-ul))*rolr*        &
                                                             UGrid%dx(i,Jsize)

        vl=v(i,Jsize);vr=v(i,Jsize+1)
        alr=(vl+vr)
        if(VCell%NEdge_ARea(i,Jsize)>=epsi) then
          rolr=((1.d0-PCell%Vof(i,Jsize)-PCell%VofS(i,Jsize))*roa/Roref+       &
                                         PCell%Vof(i,Jsize)*row/Roref)
        else
          rolr=0.d0
        end if
        if(alr>=0.d0) then
          VFlux(i,Jsize+1,2)=vl*rolr*VGrid%dx(i,Jsize)
        else
          VFlux(i,Jsize+1,2)=vr*rolr*VGrid%dx(i,Jsize)
        end if
        flux(i,Jsize+1,2)=0.5d0*(vr*vr+vl*vl-dabs(alr)*(vr-vl))*rolr*          &
                                                             VGrid%dx(i,Jsize)
      end do
    end subroutine HighOrderDensityBasedConvFluxYDir

!   in this subroutine the face area should be took into account
    subroutine DensityBasedConvectiveFlux(PGrid,UGrid,VGrid,PCell,UCell,VCell, &
                                          flux,Vflux,vb,idir,jdir)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      INTEGER(kind=it4b),INTENT(IN):: idir,jdir
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT)::flux,Vflux
      REAL(KIND=dp),INTENT(IN):: vb
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp):: uw,vs,uwn,uwp,vsn,vsp
      REAL(KIND=dp):: eta,Sx,Sy
      if(idir==1) then
        do i = 1,Isize+1
          do j = 1,Jsize
            uw = 0.5d0*(u(i,j)+u(i-1,j))
      !     for UCell both for convective velocity and scalar velocity
            if(i==ibeg) then
              Flux(i,j,1)=uw**2.d0*VFlux(i,j,1)*UGrid%dy(i,j)
            elseif(i>=ibeg+Isize-1) then
              Flux(i,j,1)=uw**2.d0*VFlux(i,j,1)*UGrid%dy(i-1,j)
            else
              eta=UCell%EtaE(i-1,j)
              uw=(1.d0-eta)*u(i-1,j)+eta*u(i,j)
      !     first order upwind
              uwp=0.5d0*(uw+dabs(uw))
              uwn=0.5d0*(uw-dabs(uw))
              Flux(i,j,1)=(uwp*u(i-1,j)+uwn*u(i,j))*UCell%AlE(i-1,j)*          &
                                                    VFlux(i,j,1)*UGrid%dy(i,j)
            end if
       !    convective velocity: u, scalar advective : v
            uw=0.5d0*(u(i-1,j+1)+u(i-1,j))
            Flux(i,j,2)=0.d0
            if(i>ibeg+Isize-1) then
              if(VCell%EEdge_Area(Isize,j)>=epsi) then
                Flux(i,j,2)=uw*0.5d0*(v(i-1,j)+v(i,j))*(Vvolf(Isize,j)*        &
                          row/Roref+(1.d0-Vvolf(Isize,j)-VCell%VofS(Isize,j))* &
                          roa/Roref)*VGrid%dy(i-1,j)
              end if
            elseif(i==ibeg) then
              if(VCell%WEdge_Area(1,j)>=epsi) then
                Flux(i,j,2)=uw*0.5d0*(v(i-1,j)+v(i,j))*(Vvolf(i,j)*row/Roref+  &
                    (1.d0-Vvolf(i,j)-VCell%VofS(i,j))*roa/Roref)*VGrid%dy(i,j)
              end if
            else
              Sy=UCell%SyN(i-1,j)
              eta=dabs(VCell%FCE(i-1,j,2)+PGrid%dy(i-1,j)/2.d0-                &
                                                   UCell%Cell_Cent(i-1,j,2))/Sy
              if(dabs(eta)>=1.d0) eta=0.5d0
              uw=(1.d0-eta)*u(i-1,j)+eta*u(i-1,j+1)
           !  second order
           !   vw=(1.d0-VCell%EtaE(i-1,j))*v(i-1,j)+VCell%EtaE(i-1,j)*v(i,j)
           !   Flux(i,j,2)=uw*vw*VCell%AlE(i-1,j)*VCell%WEdge_Area(i,j)*        &
           !                                            VGrid%dy(i,j)
           !  firt order upwind
              uwp=0.5d0*(uw+dabs(uw))
              uwn=0.5d0*(uw-dabs(uw))
              if(VCell%WEdge_Area(i,j)>=epsi) then
                Flux(i,j,2)=(uwp*v(i-1,j)*(Vvolf(i-1,j)*row/Roref+             &
                            (1.d0-Vvolf(i-1,j)-VCell%vofS(i-1,j))*roa/Roref)+  &
                             uwn*v(i,j)*(Vvolf(i,j)*row/Roref+(1.d0-Vvolf(i,j)-&
                             VCell%vofS(i,j))*roa/Roref))*VGrid%dy(i,j)
              end if
            end if
          end do
        end do
      end if

      if(jdir==1) then ! Jflux
        do i=1,Isize
          do j=1,Jsize+1
            Flux(i,j,1)=0.d0
          ! Convective velocity: v, scalar advective: u
            if(j>jbeg+Jsize-1) then
              vs=0.5d0*(v(i,j-1)+v(i+1,j-1))
              if(UCell%NEdge_Area(i,Jsize)>=epsi) then
                Flux(i,j,1)=vs*0.5d0*(u(i,j)+u(i,j-1))*(Uvolf(i,Jsize)*        &
                            row/Roref+(1.d0-Uvolf(i,Jsize)-                    &
                            UCell%vofS(i,Jsize))*roa/Roref)*UGrid%dx(i,j-1)
              end if
            elseif(j==jbeg) then
              vs=0.5d0*(v(i,j-1)+v(i+1,j-1))
              if(UCell%SEdge_Area(i,1)>=epsi) then
                Flux(i,j,1)=vs*0.5d0*(u(i,j)+u(i,j-1))*(Uvolf(i,j)*row/Roref+  &
                     (1.d0-Uvolf(i,j)-UCell%vofS(i,j))*roa/Roref)*UGrid%dx(i,j)
              end if
            else
              Sx=VCell%SxE(i,j-1)
              eta=dabs(UCell%FCN(i,j-1,1)+PGrid%dx(i,j-1)/2.d0-                &
                                                   VCell%Cell_Cent(i,j-1,1))/Sx
              if(dabs(eta)>=1.d0) eta=0.5d0
              vs=(1.d0-eta)*v(i,j-1)+eta*v(i+1,j-1)
              vsp=0.5d0*(vs+dabs(vs))
              vsn=0.5d0*(vs-dabs(vs))
              if(UCell%SEdge_Area(i,j)>=epsi) then
                Flux(i,j,1)=(vsp*u(i,j-1)*(Uvolf(i,j-1)*row/Roref+             &
                            (1.d0-Uvolf(i,j-1)-UCell%vofS(i,j-1))*roa/Roref)+  &
                             vsn*u(i,j)*(Uvolf(i,j)*row/Roref+(1.d0-Uvolf(i,j)-&
                             UCell%vofS(i,j))*roa/Roref))*UGrid%dx(i,j)
              end if
            end if
          ! Convective velocity: v, scalar advective: v
            if(j>=jbeg+Jsize-1) then
              vs=0.5d0*(v(i,j)+v(i,j-1))
              Flux(i,j,2)=vs**2.d0*VFlux(i,j,2)*VGrid%dx(i,j-1)
            elseif(j==jbeg) then
              vs=0.5d0*(v(i,j)+v(i,j-1))
              Flux(i,j,2)=vs**2.d0*VFlux(i,j,2)*VGrid%dx(i,j)
            else
              eta=VCell%EtaN(i,j-1)
              vs=(1.d0-eta)*v(i,j-1)+eta*v(i,j)
              vs=((vs-vb)*VCell%AlN(i,j-1)+vb)
              vsp=0.5d0*(vs+dabs(vs))
              vsn=0.5d0*(vs-dabs(vs))
              Flux(i,j,2)=(vsp*v(i,j-1)+vsn*v(i,j))*VFlux(i,j,2)*VGrid%dx(i,j)
            end if
          end do
        end do
      end if
    end subroutine DensityBasedConvectiveFlux

    subroutine FaceDensityFlux(PGrid,UGrid,VGrid,PCell,UCell,VCell,            &
                                      flux,Vflux,vb,dt,idir,jdir)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      INTEGER(kind=it4b),INTENT(IN):: idir,jdir
      REAL(KIND=dp),INTENT(IN):: dt
      REAL(KIND=dp),INTENT(IN):: vb
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT)::flux,VFlux
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp):: uw,vs,uwn,uwp,vsn,vsp,VofFace,diss,nxx,nyy
      REAL(KIND=dp):: eta,Sx,Sy,volf,vols
      do i = 1,Isize+idir
        do j = 1,Jsize+jdir
          if(idir==1) then
            uw=0.5d0*(u(i,j)+u(i-1,j))
            Vflux(i,j,1)=0.d0
            Flux(i,j,1)=0.d0
          ! For UCell both for convective velocity and scalar velocity
            if(i==ibeg) then
              if(UCell%WEdge_Area(1,j)>=epsi) then
                VFlux(i,j,1)=((1.d0-PCell%Vof(i,j)-PCell%VofS(i,j))*roa/Roref+ &
                                                     PCell%Vof(i,j)*row/Roref)
              end if
              Flux(i,j,1)=uw*VFlux(i,j,1)*UGrid%dy(i,j)
            elseif(i>=ibeg+Isize-1) then
              if(UCell%EEdge_Area(Isize,j)>=epsi) then
                VFlux(i,j,1)=((1.d0-PCell%Vof(Isize,j)-PCell%VofS(Isize,j))*   &
                                       roa/Roref+PCell%Vof(Isize,j)*row/Roref)
              end if
              Flux(i,j,1)=uw*VFlux(i,j,1)*UGrid%dy(i-1,j)
            else
              if(UCell%WEdge_Area(i,j)>=epsi) then
                eta=UCell%EtaE(i-1,j)
                uw=(1.d0-eta)*u(i-1,j)+eta*u(i,j)
                if(PCell%vof(i,j)>1.d0-epsi.or.PCell%vof(i,j)<epsi)then
                  VFlux(i,j,1)=((1.d0-PCell%Vof(i,j)-PCell%VofS(i,j))*roa/Roref+&
                                      PCell%Vof(i,j)*row/Roref)
                else
                  if(PCell%VofS(i,j)<epsi) then
                    call frac(PCell%nx(i,j),PCell%ny(i,j),PCell%phi(i,j)+       &
                        uw*dt/2.d0*PCell%nx(i,j),dabs(uw*dt),PGrid%dy(i,j),volf)
                    VFlux(i,j,1)=(volf*row/Roref+(1.d0-volf-vols)*roa/Roref)
                  else
                    call CellGeoCal(PCell%nxs(i,j),PCell%nys(i,j),              &
                      PCell%phis(i,j)+uw*dt/2.d0*PCell%nxs(i,j),PCell%nx(i,j),  &
                      PCell%ny(i,j),PCell%phi(i,j)+uw*dt/2.d0*PCell%nx(i,j),    &
                                    dabs(uw*dt),PGrid%dy(i,j),vols,volf)
                    VFlux(i,j,1)=(volf*row/Roref+(1.d0-volf-vols)*roa/Roref)
                  end if
                end if
              end if
              Flux(i,j,1)=uw*UCell%AlE(i-1,j)*VFlux(i,j,1)*UGrid%dy(i,j)
            end if
          ! Convective velocity: u, scalar advective : v
            uw = 0.5d0*(u(i-1,j+1)+u(i-1,j))
            VFlux(i,j,2)=0.d0
            Flux(i,j,2)=0.d0
            if(i>ibeg+Isize-1) then
              if(VCell%EEdge_Area(Isize,j)>=epsi) then
                Flux(i,j,2)=uw*(Vvolf(Isize,j)*row/Roref+                      &
                    (1.d0-Vvolf(Isize,j)-VCell%Vofs(Isize,j))*roa/Roref)*      &
                                                             VGrid%dy(i-1,j)
              end if
            elseif(i==ibeg) then
              if(VCell%WEdge_ARea(1,j)>=epsi) then
                Flux(i,j,2)=uw*(Vvolf(i,j)*row/Roref+                          &
                    (1.d0-Vvolf(i,j)-VCell%vofS(i,j))*roa/Roref)*VGrid%dy(i,j)
              end if
            else
              if(VCell%WEdge_Area(i,j)>=epsi) then
                Sy=UCell%SyN(i-1,j)
                eta=dabs(VCell%FCE(i-1,j,2)+PGrid%dy(i-1,j)/2.d0-              &
                                                UCell%Cell_Cent(i-1,j,2))/Sy
                if(dabs(eta)>=1.d0) eta=0.5d0
                uw=(1.d0-eta)*u(i-1,j)+eta*u(i-1,j+1)
                uwp=0.5d0*(uw+dabs(uw))
                uwn=0.5d0*(uw-dabs(uw))
                Flux(i,j,2)=(uwp*(Vvolf(i-1,j)*row/Roref+(1.d0-Vvolf(i-1,j)-   &
                             VCell%vofS(i-1,j))*roa/Roref)+uwn*(Vvolf(i,j)*    &
                             row/Roref+(1.d0-Vvolf(i,j)-VCell%vofS(i,j))*      &
                             roa/Roref))*VGrid%dy(i,j)
              end if
            end if
          end if
          if(jdir==1) then ! Jflux
            VFlux(i,j,1)=0.d0
            Flux(i,j,1)=0.d0
        ! Convective velocity: v, scalar advective: u
            vs = 0.5d0*(v(i,j-1)+v(i+1,j-1))
            if(j>jbeg+Jsize-1) then
              if(UCell%NEdge_Area(i,Jsize)>=epsi) then
                Flux(i,j,1)=vs*(Uvolf(i,Jsize)*row/Roref+                      &
                        (1.d0-Uvolf(i,Jsize)-UCell%vofS(i,Jsize))*roa/Roref)*  &
                                                             UGrid%dx(i,j-1)
              end if
            elseif(j==jbeg) then
              if(UCell%SEdge_Area(i,j)>=epsi) then
                Flux(i,j,1)=vs*(Uvolf(i,j)*row/Roref+                          &
                  (1.d0-Uvolf(i,j)-UCell%vofS(i,j))*roa/Roref)*UGrid%dx(i,j)
              end if
            else
              if(UCell%SEdge_Area(i,j)>=epsi) then
                Sx=VCell%SxE(i,j-1)
                eta=dabs(UCell%FCN(i,j-1,1)+PGrid%dx(i,j-1)/2.d0-              &
                                            VCell%Cell_Cent(i,j-1,1))/Sx
                if(dabs(eta)>=1.d0) eta=0.5d0
                vs=(1.d0-eta)*v(i,j-1)+eta*v(i+1,j-1)
                vsp = 0.5d0*(vs+dabs(vs))
                vsn = 0.5d0*(vs-dabs(vs))
                Flux(i,j,1)=(vsp*(Uvolf(i,j-1)*row/Roref+(1.d0-Uvolf(i,j-1)-   &
                             UCell%vofS(i,j-1))*roa/Roref)+vsn*(Uvolf(i,j)*    &
                             row/Roref+(1.d0-Uvolf(i,j)-UCell%vofS(i,j))*      &
                             roa/Roref))*UGrid%dx(i,j)
              end if
            end if
            vs = 0.5d0*(v(i,j)+v(i,j-1))
            Vflux(i,j,2)=0.d0
            Flux(i,j,2)=0.d0
            if(j>=jbeg+Jsize-1) then
              if(VCell%NEdge_Area(i,Jsize)>=epsi) then
                VFlux(i,j,2)=((1.d0-PCell%Vof(i,Jsize)-PCell%vofS(i,Jsize))*   &
                               roa/Roref+PCell%Vof(i,Jsize)*row/Roref)
              end if
              Flux(i,j,2)=vs*VFlux(i,j,2)*VGrid%dx(i,j-1)
            elseif(j==jbeg) then
              if(VCell%SEdge_Area(i,j)>=epsi) then
                VFlux(i,j,2)=((1.d0-PCell%Vof(i,j)-PCell%vofS(i,j))*roa/Roref+ &
                                    PCell%Vof(i,j)*row/Roref)
              end if
              Flux(i,j,2)=vs*VFlux(i,j,2)*VGrid%dx(i,j)
            else
              if(VCell%SEdge_Area(i,j)>=epsi) then
                eta=VCell%EtaN(i,j-1)
                vs=(1.d0-eta)*v(i,j-1)+eta*v(i,j)
                vs=((vs-vb)*VCell%AlN(i,j-1)+vb)
                if(PCell%vof(i,j)>1.d0-epsi.or.PCell%vof(i,j)<epsi)then
                  VFlux(i,j,2)=(PCell%vof(i,j)*row/Roref+                      &
                               (1.d0-PCell%vof(i,j)-PCell%vofS(i,j))*roa/Roref)
                else
                  if(PCell%vofS(i,j)<epsi) then
                    call frac(PCell%nx(i,j),PCell%ny(i,j),PCell%phi(i,j)+      &
                        vs*dt/2.d0*PCell%ny(i,j),PGrid%dx(i,j),dabs(vs*dt),volf)
                    VFlux(i,j,2)=(volf*row/Roref+(1.d0-volf-vols)*roa/Roref)
                  else
                    call CellGeoCal(PCell%nxs(i,j),PCell%nys(i,j),             &
                       PCell%phis(i,j)+vs*dt/2.d0*PCell%nys(i,j),PCell%nx(i,j),&
                       PCell%ny(i,j),PCell%phi(i,j)+vs*dt/2.d0*PCell%ny(i,j),  &
                                     PGrid%dx(i,j),dabs(vs*dt),vols,volf)
                    VFlux(i,j,2)=(volf*row/Roref+(1.d0-volf-vols)*roa/Roref)
                  end if
                end if
                Flux(i,j,2)=vs*Vflux(i,j,2)*VGrid%dx(i,j)
              end if
            end if
          end if
        end do
      end do
    end subroutine FaceDensityFlux

    SUBROUTINE DiffusiveFlux(PGrid,UGrid,VGrid,PCell,UCell,VCell,flux,Eflux,   &
                                                     ExEFlux,vb,idir,jdir)
      IMPLICIT NONE
      INTEGER(kind=it4b),INTENT(IN):: idir,jdir
      REAL(KIND=dp),INTENT(IN):: vb
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT)::flux,Eflux,ExEFlux
      REAL(KIND=dp):: tol,VofFace,VofSFace,Sx,Sy
      tol = 1.d-24
      Eflux(:,:,:)=0.d0
      EXEFlux(:,:,:)=0.d0
      if(idir==1) then
        do j=1,Jsize
          flux(1,j,1)=((1.d0-Uvolf(1,j)/(1.d0-UCell%VofS(1,j)+tol))*nua/nuref+ &
             Uvolf(1,j)/(1.d0-UCell%VofS(1,j)+tol)*nuw/nuref)*                 &
             UCell%WEdge_Area(1,j)*UGrid%dy(1,j)/PGrid%dx(1,j)/Rey
          flux(1,j,2)=((1.d0-Vvolf(1,j)/(1.d0-VCell%VofS(1,j)+tol))*nua/nuref+ &
             Vvolf(1,j)/(1.d0-VCell%vofS(1,j)+tol)*nuw/nuref)*                 &
             VCell%WEdge_Area(1,j)*VGrid%dy(1,j)/PGrid%dx(1,j)/Rey
          do i=2,Isize
            Sx=UCell%SxE(i-1,j)
            Sy=UCell%Cell_Cent(i,j,2)-UCell%Cell_Cent(i-1,j,2)
            VofFace=0.5d0*(Uvolf(i,j)/(1.d0-UCell%VofS(i,j)+tol)+              &
                           Uvolf(i-1,j)/(1.d0-UCell%VofS(i-1,j)+tol))
            flux(i,j,1)=((1.d0-VofFace)*nua/nuref+VofFace*nuw/nuref)/Rey*      &
                         UCell%WEdge_Area(i,j)*UGrid%dy(i,j)/Sx
            if(dabs(Sy)>=1.d-3*UGrid%dy(i,j)) then
              EFlux(i,j,1)=((1.d0-VofFace)*nua/nuref+VofFace*nuw/nuref)/Rey*   &
                         UCell%WEdge_Area(i,j)*UCell%DAlE(i-1,j)*UGrid%dy(i,j)
            end if
            Sx=VCell%SxE(i-1,j)
            Sy=VCell%Cell_Cent(i,j,2)-VCell%Cell_Cent(i-1,j,2)
            VofFace=0.5d0*(Vvolf(i,j)/(1.d0-VCell%VofS(i,j)+tol)+              &
                           Vvolf(i-1,j)/(1.d0-VCell%VofS(i-1,j)+tol))
            flux(i,j,2)=((1.d0-VofFace)*nua/nuref+VofFace*nuw/nuref)/Rey*      &
                         VCell%WEdge_Area(i,j)*VGrid%dy(i,j)/Sx
            if(dabs(Sy)>=1.d-3*VGrid%dy(i,j)) then
              EFlux(i,j,2)=((1.d0-VofFace)*nua/nuref+VofFace*nuw/nuref)/Rey*   &
                         VCell%WEdge_Area(i,j)*VCell%DAlE(i-1,j)*VGrid%dy(i,j)
              ExEFlux(i,j,2)=vb*EFlux(i,j,2)
            end if
          end do
          flux(Isize+1,j,1)=((1.d0-Uvolf(Isize,j)/                             &
                 (1.d0-UCell%VofS(Isize,j)+tol))*nua/nuref+Uvolf(Isize,j)/     &
                 (1.d0-UCell%VofS(Isize,j)+tol)*nuw/nuref)/Rey*                &
                 UCell%EEdge_Area(Isize,j)*UGrid%dy(Isize,j)/PGrid%dx(Isize,j)
          flux(Isize+1,j,2)=((1.d0-Vvolf(Isize,j)/                             &
                 (1.d0-VCell%VofS(Isize,j)+tol))*nua/nuref+Vvolf(Isize,j)/     &
                 (1.d0-VCell%vofS(Isize,j)+tol)*nuw/nuref)/Rey*                &
                 VCell%EEdge_Area(Isize,j)*VGrid%dy(Isize,j)/PGrid%dx(Isize,j)
        end do
      elseif(jdir==1) then
        do i=1,Isize
          flux(i,1,1)=((1.d0-Uvolf(i,1)/(1.d0-UCell%VofS(i,1)+tol))*nua/nuref+ &
             Uvolf(i,1)/(1.d0-UCell%VofS(i,1)+tol)*nuw/nuref)/Rey*             &
             UCell%SEdge_Area(i,1)*UGrid%dx(i,1)/PGrid%dy(i,1)
          flux(i,1,2)=((1.d0-Vvolf(i,1)/(1.d0-VCell%VofS(i,1)+tol))*nua/nuref+ &
             Vvolf(i,1)/(1.d0-VCell%vofS(i,1)+tol)*nuw/nuref)/Rey*             &
             VCell%SEdge_Area(i,1)*VGrid%dx(i,1)/PGrid%dy(i,1)
          do j=2,Jsize
            Sx=UCell%Cell_Cent(i,j,1)-UCell%Cell_Cent(i,j-1,1)
            Sy=UCell%SyN(i,j-1)
            VofFace=0.5d0*(Uvolf(i,j-1)/(1.d0-UCell%VofS(i,j-1)+tol)+          &
                           Uvolf(i,j)/(1.d0-UCell%VofS(i,j)+tol))
            flux(i,j,1)=((1.d0-VofFace)*nua/nuref+VofFace*nuw/nuref)/Rey*      &
                         UCell%SEdge_Area(i,j)*UGrid%dx(i,j)/Sy
            if(dabs(Sx)>=1.d-3*UGrid%dx(i,j)) then
              EFlux(i,j,1)=((1.d0-VofFace)*nua/nuref+VofFace*nuw/nuref)/Rey*   &
                         UCell%SEdge_Area(i,j)*UCell%DAlN(i,j-1)*UGrid%dx(i,j)
            end if
            Sx=VCell%Cell_Cent(i,j,1)-VCell%Cell_Cent(i,j-1,1)
            Sy=VCell%SyN(i,j-1)
            VofFace=0.5d0*(Vvolf(i,j-1)/(1.d0-VCell%VofS(i,j-1)+tol)+          &
                           Vvolf(i,j)/(1.d0-VCell%VofS(i,j)+tol))
            flux(i,j,2)=((1.d0-VofFace)*nua/nuref+VofFace*nuw/nuref)/Rey*      &
                         VCell%SEdge_Area(i,j)*VGrid%dx(i,j)/Sy
            if(dabs(Sx)>=1.d-3*VGrid%dx(i,j)) then
              EFlux(i,j,2)=((1.d0-VofFace)*nua/nuref+VofFace*nuw/nuref)/Rey*   &
                         VCell%SEdge_Area(i,j)*VCell%DAlN(i,j-1)*VGrid%dx(i,j)
              ExEFlux(i,j,2)=vb*EFlux(i,j,2)
            end if
          end do
          flux(i,Jsize+1,1)=((1.d0-Uvolf(i,Jsize)/                             &
                   (1.d0-UCell%VofS(i,Jsize)+tol))*nua/nuref+Uvolf(i,Jsize)/   &
                   (1.d0-UCell%VofS(i,Jsize)+tol)*nuw/nuref)/Rey*              &
                   UCell%NEdge_Area(i,Jsize)*UGrid%dx(i,Jsize)/PGrid%dy(i,Jsize)
          flux(i,Jsize+1,2)=((1.d0-Vvolf(i,Jsize)/                             &
                   (1.d0-VCell%VofS(i,Jsize)+tol))*nua/nuref+Vvolf(i,Jsize)/   &
                   (1.d0-VCell%vofS(i,Jsize)+tol)*nuw/nuref)/Rey*              &
                   VCell%NEdge_Area(i,Jsize)*VGrid%dx(i,Jsize)/PGrid%dy(i,Jsize)
        end do
      end if
    end subroutine DiffusiveFlux

    subroutine Mixing_Procedure(TGrid,TCell,idir,jdir,flux_mix)
      TYPE(Grid),INTENT(IN):: TGrid
      TYPE(Cell),INTENT(IN):: TCell
      INTEGER(kind=it4b),INTENT(IN):: idir,jdir
      REAL(KIND=dp),DIMENSION(:,:),INTENT(INOUT):: flux_mix
      INTEGER(kind=it4b):: i,j,ii,jj
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: flux
      REAL(KIND=dp):: Fmix(3),Beta(3),Vtgt(3),temp
      allocate(flux(ibeg:Isize,jbeg:Jsize))
      flux(:,:)=flux_mix(:,:)
      do i=ibeg,Isize-idir
        do j=jbeg,Jsize-jdir
          if(TCell%VofS(i,j)>=0.5d0+epsi.and.TCell%Posnu(i,j)/=-1) then
            ii=i+INT(sign(1.d0,TCell%nxS(i,j)));jj=j
            if(ii<=Isize-idir.and.ii>=1) then
              Vtgt(1)=1.d0-TCell%VofS(ii,jj)
              Fmix(1)=(1.d0-TCell%vofS(i,j))*Vtgt(1)*(flux(ii,jj)-flux(i,j))
              Beta(1)=TCell%nxS(i,j)**2.d0*Vtgt(1)**5.d0
            else
              Fmix(1)=0.d0;Beta(1)=0.d0;Vtgt(1)=1.d0
            end if
            ii=i;jj=j+INT(sign(1.d0,TCell%nyS(i,j)))
            if(jj<=Jsize-jdir.and.jj>=1) then
              Vtgt(2)=1.d0-TCell%VofS(ii,jj)
              Fmix(2)=(1.d0-TCell%VofS(i,j))*Vtgt(2)*(flux(ii,jj)-flux(i,j))
              Beta(2)=TCell%nyS(i,j)**2.d0*Vtgt(2)**5.d0
            else
              Fmix(2)=0.d0;Beta(2)=0.d0;Vtgt(2)=1.d0
            end if
            ii=i+INT(sign(1.d0,TCell%nxS(i,j)))
            jj=j+INT(sign(1.d0,TCell%nyS(i,j)))
            if(ii<=Isize-idir.and.ii>=1.and.jj<=Jsize-jdir.and.jj>=1) then
              Vtgt(3)=1.d0-TCell%VofS(ii,jj)
              Fmix(3)=(1.d0-TCell%vofS(i,j))*Vtgt(3)*(flux(ii,jj)-flux(i,j))
              Beta(3)=dabs(TCell%nxS(i,j)*TCell%nyS(i,j))*Vtgt(3)**5.d0
            else
              Fmix(3)=0.d0;Beta(3)=0.d0;Vtgt(3)=1.d0
            end if
            temp=sum(Beta)
            do ii=1,3
              Beta(ii)=Beta(ii)/temp
              Fmix(ii)=Fmix(ii)*Beta(ii)/(Beta(ii)*(1.d0-TCell%VofS(i,j))+     &
                                                                   Vtgt(ii))
            end do
            flux_mix(i,j)=flux(i,j)+1.d0/(1.d0-TCell%VofS(i,j))*(sum(Fmix))
            ii=i+INT(sign(1.d0,TCell%nxS(i,j)));jj=j
            if(ii<=Isize-idir.and.ii>=1) then    ! for target cell
              flux_mix(ii,jj)=flux_mix(ii,jj)-1.d0/Vtgt(1)*Fmix(1)
            end if
            ii=i;jj=j+INT(sign(1.d0,TCell%nyS(i,j)))
            if(jj<=Jsize-jdir.and.jj>=1) then    ! for target cell
              flux_mix(ii,jj)=flux_mix(ii,jj)-1.d0/Vtgt(2)*Fmix(2)
            end if
            ii=i+INT(sign(1.d0,TCell%nxS(i,j)))
            jj=j+INT(sign(1.d0,TCell%nyS(i,j)))
            if(ii<=Isize-idir.and.ii>=1.and.jj<=Jsize-jdir.and.jj>=1) then ! for target cell
              flux_mix(ii,jj)=flux_mix(ii,jj)-1.d0/Vtgt(3)*Fmix(3)
            end if
            if(isnan(Flux_Mix(i,j)).or.dabs(Flux_Mix(i,j))>1.d4.or.            &
               isnan(Flux_Mix(ii,jj)).or.dabs(Flux_Mix(ii,jj))>1.d4) then
              print*,i,j
              print*,ii,jj
              print*,Fmix(1),Fmix(2),Fmix(3)
              print*,Flux_Mix(i,j),Flux_Mix(ii,jj)
              print*,'MixingProcedure_PredictorUV_2169'
            end if
          end if
        end do
      end do
      deallocate(flux)
    end subroutine Mixing_Procedure

    subroutine PredictorVelocityBoundaryCondition(Pred,TVar)
      IMPLICIT NONE
      TYPE(Predictor),INTENT(INOUT):: Pred
      TYPE(Variables),INTENT(IN):: TVar
      REAL(KIND=dp),PARAMETER:: Twall = 300.d0
      INTEGER(kind=it4b):: i,j
      do j = jbeg,Jsize+jbeg-1
      ! Outlet
        Pred%u(Isize+ibeg-1,j)=Pred%u(Isize+ibeg-2,j)
      end do
      do i = ibeg,Isize+ibeg-1
      ! Open air
        Pred%v(i,Jsize+jbeg-1)=Pred%v(i,Jsize+jbeg-2)
      end do
    end subroutine PredictorVelocityBoundaryCondition

    SUBROUTINE frac(nx,ny,diss,dx,dy,vrt)
       IMPLICIT NONE
       REAL(KIND=dp):: slop_eps,nx,ny,diss,vrt,tnx,tny,dx,dy
       REAL(KIND=dp):: xx,yy,topvf,rightvf,totalarea
       slop_eps = 1.d-20
       tnx = dabs(nx)
       tny = dabs(ny)
       vrt = 0.d0
       if(diss+0.5d0*(tnx*dx+tny*dy)<=0.d0) then
          vrt = 1.d0
          return
       end if
       if(diss-0.5d0*(tnx*dx+tny*dy)>=0.d0) then
          vrt = 0.d0
          return
       end if
       if(tnx <= slop_eps) then
          vrt = 0.5d0-diss/dy
          return
       end if
       if(tny <= slop_eps) then
          vrt = 0.5d0-diss/dx
          return
       end if
       xx = (0.5d0*tny*dy-diss)/dx/tnx+0.5d0
       yy = (0.5d0*tnx*dx-diss)/dy/tny+0.5d0
       totalarea = 0.5d0*xx*yy
       topvf = 0.d0
       if(yy > 1.d0) topvf =((yy-1.d0)/yy)**2.d0
       rightvf = 0.d0
       if(xx > 1.d0) rightvf = ((xx-1.d0)/xx)**2.d0
       vrt = totalarea*(1.d0-topvf-rightvf)
       if(vrt<0.d0) then
         print*,vrt,diss,nx,ny
         print*,tnx,tny,xx,yy
         print*,dy,dx
         write(*,*) 'vrt less than 0'
         print*,topvf,rightvf
         pause 'PredictorUV_1068'
       end if
       return
    END SUBROUTINE

    SUBROUTINE CellGeoCal(nxs,nys,phis,nxl,nyl,phil,dx,dy,vols,volf)
        IMPLICIT NONE
        REAL(KIND=dp),INTENT(IN):: nxs,nys,phis,nxl,nyl,phil,dx,dy
        REAL(KIND=dp),INTENT(OUT):: vols,volf
        REAL(KIND=dp),DIMENSION(:,:):: node(6,2),nodel(7,2)
        INTEGER(kind=it4b):: temp,templ,k
        REAL(KIND=dp):: dx2,dy2,pdt(6),dpt(4),vol,epsil
        epsil=1.d-24
        dx2 = dx/2.d0
        dy2 = dy/2.d0
        temp = 1
        node = 0.d0
     !  calculate the distance of all cell's nodes to interface
        dpt(1) =  dx2*nxs-dy2*nys+phis
        dpt(2) = -dx2*nxs-dy2*nys+phis
        dpt(3) = -dx2*nxs+dy2*nys+phis
        dpt(4) =  dx2*nxs+dy2*nys+phis
     !  find nodes those belong to fluid field(including the cutting nodes between
     !  interface and cell edges)
        If(dpt(1)>=0.d0) then
          node(temp,1) = dx2
          node(temp,2) = -dy2
          temp = temp+1
        End if
        If(dpt(1)*dpt(2)<0.d0) then
          node(temp,1) = (dy2*nys-phis)/nxs
          node(temp,2) = -dy2
          temp = temp+1
        End if
        If(dpt(2)>=0.d0) then
          node(temp,1) = -dx2
          node(temp,2) = -dy2
          temp = temp+1
        End if
        If(dpt(2)*dpt(3)<0.d0) then
          node(temp,1) = -dx2
          node(temp,2) = (dx2*nxs-phis)/nys
          temp = temp+1
        End if
        If(dpt(3)>=0.d0) then
          node(temp,1) = -dx2
          node(temp,2) = dy2
          temp = temp+1
        End if
        If(dpt(3)*dpt(4)<0.d0) then
          node(temp,1) = (-dy2*nys-phis)/nxs
          node(temp,2) = dy2
          temp = temp+1
        End if
        If(dpt(4)>=0.d0) then
          node(temp,1) = dx2
          node(temp,2) = dy2
          temp = temp+1
        End if
        If(dpt(4)*dpt(1)<0.d0) then
          node(temp,1) = dx2
          node(temp,2) = (-dx2*nxs-phis)/nys
          temp = temp+1
        End if
        node(temp,1) = node(1,1)
        node(temp,2) = node(1,2)
        vol = 0.d0
    !   applying Gauss theorem to find the volume of both fluid and gas
        do k = 1,temp-1
          vol = vol+0.5d0*(node(k,1)*node(k+1,2)-node(k+1,1)*node(k,2))
          pdt(k)=node(k,1)*nxl+node(k,2)*nyl+phil
        end do
    !   Volume fraction of solid in given cell
        vols = 1.d0-dabs(vol/(dx*dy))
        pdt(temp)=pdt(1)
    !   find the volume of liquid inside cell
        templ=1
        do k=1,temp-1
          if(pdt(k)<=0.d0) then
            nodel(templ,1)=node(k,1)
            nodel(templ,2)=node(k,2)
            templ=templ+1
          end if
          if(pdt(k)*pdt(k+1)<0.d0) then
            nodel(templ,1)=node(k,1)+(node(k+1,1)-node(k,1))*dabs(pdt(k))/      &
                                                   (dabs(pdt(k))+dabs(pdt(k+1)))
            nodel(templ,2)=node(k,2)+(node(k+1,2)-node(k,2))*dabs(pdt(k))/      &
                                                   (dabs(pdt(k))+dabs(pdt(k+1)))
            templ=templ+1
          end if
        end do

        nodel(templ,1) = nodel(1,1)
        nodel(templ,2) = nodel(1,2)
        vol = 0.d0
        do k = 1,templ-1
          vol = vol+0.5d0*(nodel(k,1)*nodel(k+1,2)-nodel(k+1,1)*nodel(k,2))
        end do
        volf=dabs(vol/(dx*dy))
    end subroutine CellGeoCal

    function MUSCLLimiter(x,opt,tolim) result (y)
      REAL(KIND=dp),INTENT(IN):: x
      INTEGER(kind=it4b),INTENT(IN):: opt
      REAL(KIND=dp),INTENT(IN):: tolim
      REAL(KIND=dp):: y
      select case(opt)
      ! MinMod limiter
        case (1)
          y=dmax1(0.d0,dmin1(tolim,x))
      ! SuperBee limiter
        case (2)
          y=dmax1(0.d0,dmin1(tolim*x,1.d0),dmin1(x,tolim))
      end select
    end function MUSCLLimiter
end module PredictorUV
