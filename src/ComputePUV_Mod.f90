Module ComputePUV
 !! Description:
 !! The module calculates the pressure and velocities (correction step). lied.
 ! Current Code Owner: SIMCOFlow
 ! Code Description:
 ! Language: Fortran 90.
 ! Software Standards: "European Standards for Writing and
 ! Documenting Exchangeable Fortran 90 Code".
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! Author : Son Tung Dang
 !        : NTNU,SINTEF
 ! Date : 20.09.2019  
    USE PrecisionVar
    USE Mesh, ONLY : Grid, Cell,ibeg, jbeg, Isize, Jsize, ight, jght
    USE Clsvof, ONLY : SolidObject
    USE StateVariables, ONLY : TVariables
    USE Constants, ONLY : epsi, roa, row
    USE PredictorUV, ONLY : Predictor, PoissonCoefficient, Predictor_UV
    USE ProjectionP, ONLY : Projection, PoissonEquationSolver
    USE Particles, ONLY : TParticle
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: UpdatePUV
    Interface UpdatePUV
      Module Procedure UpdatePUV
    End interface
    Contains
    Subroutine UpdatePUV(UGrid,VGrid,PGrid,PCellO,UCellO,VCellO,PCell,UCell,    &
           VCell,TVar,Flux_n1,TraPar,BoomCase,dt,itt)
      !! The subroutine computes the flow field parameters in the next time step.
      !!   
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN)   :: UGrid,VGrid,PGrid
        !! The input grid
      TYPE(Cell),INTENT(INOUT):: PCell,UCell,VCell
        !! The input and output cell
      TYPE(Cell),INTENT(INOUT):: PCellO,UCellO,VCellO
        !! The input and output cell a previous time step
      TYPE(Variables),INTENT(INOUT):: TVar
        !! The state variables
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT):: Flux_n1
        !! The total flux at the previous time step
      TYPE(Particle),INTENT(INOUT):: TraPar
        !! The particle information
      TYPE(SolidObject),INTENT(IN):: BoomCase
        !! The boom
      REAL(KIND=dp),INTENT(IN):: dt
        !! The time step 
      INTEGER(kind=it8b),INTENT(IN):: itt
      TYPE(PoissonCoefficient):: PU,PV
      TYPE(Predictor):: Pred
      TYPE(Projection):: Proj
      INTEGER(kind=it4b):: i,j,ii,jj
      REAL(KIND=dp):: dps
      REAL(KIND=dp):: BetaP,BetaM,BetaW,BetaD,Yint,Hjump,Lamda,tol
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: GradPUV
      allocate(Pred%u(ibeg-ight:ibeg+ight+Isize-1,jbeg-jght:jbeg+jght+Jsize-1))
      allocate(Pred%v(ibeg-ight:ibeg+ight+Isize-1,jbeg-jght:jbeg+jght+Jsize-1))
      allocate(Proj%Pp(ibeg-ight:ibeg+ight+Isize-1,jbeg-jght:jbeg+jght+Jsize-1))
      allocate(PU%Dp(ibeg-ight:jbeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
      allocate(PV%Dp(ibeg-ight:jbeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
      allocate(GradPUV(Isize,Jsize))
      BetaP = 1.d0/(row/TVar%Roref)
      BetaM = 1.d0/(roa/TVar%Roref)
      Pred%u(:,:) = TVar%u(:,:)
      Pred%v(:,:) = TVar%v(:,:)
      tol=1.d-24

      Proj%Pp(:,:) = 0.d0 
      dps=0.d0
      call Predictor_UV(PGrid,UGrid,VGrid,PCellO,UCellO,VCellO,PCell,UCell,    &
                              VCell,TVar,PU,PV,Pred,Flux_n1,TraPar,VolParU,    &
                              BoomCase,dt,itt)

      call PoissonEquationSolver(PGrid,UGrid,VGrid,PCellO,PCell,UCell,VCell,   &
                                 TVar,Pred,PU,PV,Proj,BoomCase%vs,dt,itt)
      ! Correct u
      do i=1,Isize-1
        do j=1,Jsize
          GRadPUV(i,j)=0.d0
          if(UCell%Posnu(i,j)/=-1.and.UCell%MoExCell(i,j)/=1) then
            if((PCell%vof(i,j)>=0.5d0.and.PCell%vofS(i,j)<epsi).or.            &
                (PCell%phi(i,j)<0.d0.and.PCell%vofS(i,j)>=epsi))then 
              if((PCell%vof(i+1,j)<0.5d0.and.PCell%vofS(i+1,j)<epsi).or.       &
                 (PCell%phi(i+1,j)>=0.d0.and.PCell%vofS(i+1,j)>=epsi))then
                Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+              &
                                            dabs(PCell%phi(i+1,j))+tol)
                BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP
                GradPUV(i,j)=(Proj%Pp(i+1,j)-Proj%Pp(i,j))*BetaM*BetaP/BetaW*  &
                              PU%Dp(i,j)/UGrid%dx(i,j)
                TVar%u(i,j)=Pred%u(i,j)-GradPUV(i,j)
              elseif((PCell%vof(i+1,j)>=0.5d0.and.PCell%vofS(i+1,j)<epsi).or.  &
                     (PCell%phi(i+1,j)<0.d0.and.PCell%vofS(i+1,j)>=epsi))then
                BetaW=BetaM
                GradPUV(i,j)=(Proj%Pp(i+1,j)-Proj%Pp(i,j))*BetaM*BetaP/BetaW*  &
                              PU%Dp(i,j)/UGrid%dx(i,j)
                TVar%u(i,j)=Pred%u(i,j)-GradPUV(i,j)
              end if
            elseif((PCell%vof(i,j)<0.5d0.and.PCell%vofS(i,j)<epsi).or.         &
                   (PCell%phi(i,j)>=0.d0.and.PCell%vofS(i,j)>=epsi))then
              if((PCell%vof(i+1,j)>=0.5d0.and.PCell%vofS(i+1,j)<epsi).or.      &
                 (PCell%phi(i+1,j)<0.d0.and.PCell%vofS(i+1,j)>=epsi)) then
                Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+              &
                                          dabs(PCell%phi(i+1,j))+tol)
                BetaW=Lamda*BetaP+(1.d0-Lamda)*BetaM
                GradPUV(i,j)=(Proj%Pp(i+1,j)-Proj%Pp(i,j))*BetaM*BetaP/BetaW*  &
                              PU%Dp(i,j)/UGrid%dx(i,j)
                TVar%u(i,j)=Pred%u(i,j)-GradPUV(i,j)
              elseif((PCell%vof(i+1,j)<0.5d0.and.PCell%vofS(i+1,j)<epsi).or.   &
                     (PCell%phi(i+1,j)>=0.d0.and.PCell%vofS(i+1,j)>=epsi))then
                BetaW=BetaP
                GradPUV(i,j)=(Proj%Pp(i+1,j)-Proj%Pp(i,j))*BetaM*BetaP/BetaW*  &
                              PU%Dp(i,j)/UGrid%dx(i,j)
                TVar%u(i,j)=Pred%u(i,j)-GradPUV(i,j)
              end if
            end if
          end if
        end do
      end do

    ! Cell-Linking method for small cell
      do i=1,Isize-1
        do j=1,Jsize
          if(UCell%MoExCell(i,j)==1) then
            ii=UCell%MsCe(i,j,1)
            jj=UCell%MsCe(i,j,2)
       !     TVar%u(i,j)=0.d0 !Pred%u(i,j)!-GradPUV(ii,jj)/PU%Dp(ii,jj)*PU%Dp(i,j)
          end if
        end do
      end do
      ! Correct v
      do i=1,Isize
        do j=1,Jsize-1
          GRadPUV(i,j)=0.d0
          if(VCell%Posnu(i,j)/=-1.and.VCell%MoExCell(i,j)/=1) then
            if((PCell%vof(i,j)>=0.5d0.and.PCell%vofS(i,j)<epsi).or.            &
                (PCell%phi(i,j)<0.d0.and.PCell%vofS(i,j)>=epsi))then  ! this cell is wet
              if((PCell%vof(i,j+1)<0.5d0.and.PCell%vofS(i,j+1)<epsi).or.       &
                 (PCell%phi(i,j+1)>=0.d0.and.PCell%vofS(i,j+1)>=epsi))then
                Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+              &
                                            dabs(PCell%phi(i,j+1))+tol)
                BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP
                GradPUV(i,j)=(Proj%Pp(i,j+1)-Proj%Pp(i,j))*BetaM*BetaP/BetaW*  &
                              PV%Dp(i,j)/VGrid%dy(i,j)
                TVar%v(i,j)=Pred%v(i,j)-GradPUV(i,j)
              elseif((PCell%vof(i,j+1)>=0.5d0.and.PCell%vofS(i,j+1)<epsi).or.  &
                     (PCell%phi(i,j+1)<0.d0.and.PCell%vofS(i,j+1)>=epsi))then
                BetaW=BetaM
                GRadPUV(i,j)=(Proj%Pp(i,j+1)-Proj%Pp(i,j))*BetaM*BetaP/BetaW*  &
                              PV%Dp(i,j)/VGrid%dy(i,j)
                TVar%v(i,j)=Pred%v(i,j)-GradPUV(i,j)
              end if
            elseif((PCell%vof(i,j)<0.5d0.and.PCell%vofS(i,j)<epsi).or.         &
                   (PCell%phi(i,j)>=0.d0.and.PCell%vofS(i,j)>=epsi))then
              if((PCell%vof(i,j+1)>=0.5d0.and.PCell%vofS(i,j+1)<epsi).or.      &
                 (PCell%phi(i,j+1)<0.d0.and.PCell%vofS(i,j+1)>=epsi))then
                Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+              &
                                            dabs(PCell%phi(i,j+1))+tol)
                BetaW=Lamda*BetaP+(1.d0-Lamda)*BetaM
                GradPUV(i,j)=(Proj%Pp(i,j+1)-Proj%Pp(i,j))*BetaM*BetaP/BetaW*  &
                              PV%Dp(i,j)/VGrid%dy(i,j)
                TVar%v(i,j)=Pred%v(i,j)-GradPUV(i,j)
              elseif((PCell%vof(i,j+1)<0.5d0.and.PCell%vofS(i,j+1)<epsi).or.   &
                 (PCell%phi(i,j+1)>=0.d0.and.PCell%vofS(i,j+1)>=epsi))then
                BetaW=BetaP
                GradPUV(i,j)=(Proj%Pp(i,j+1)-Proj%Pp(i,j))*BetaM*BetaP/BetaW*  &
                              PV%Dp(i,j)/VGrid%dy(i,j)
                TVar%v(i,j)=Pred%v(i,j)-GradPUV(i,j)
              end if
            end if
          end if
        end do
      end do

      ! Cell-Linking method for small cell
      do i=1,Isize
        do j=1,Jsize-1
          if(VCell%MoExCell(i,j)==1) then
            ii=VCell%MsCe(i,j,1)
            jj=VCell%MsCe(i,j,2)
          end if
        end do
      end do
      do i = 1,Isize
        do j = 1,Jsize
          if(PCell%Posnu(i,j)/=-1) then
            TVar%p(i,j) = Proj%Pp(i,j)!+TVar%p(i,j)
          end if
        end do
      end do
!      do j = 1,Jsize
!        TVar%u(Isize,j)=((-TVar%v(Isize,j)*PCell%NEdge_Area(Isize,j)+          &
!                          TVar%v(Isize,j-1)*PCell%SEdge_Area(Isize,j))*        &
!                          PGrid%dx(Isize,j)+TVar%u(Isize-1,j)*                 &
!                          PCell%WEdge_Area(Isize,j)*PGrid%dy(Isize,j))/        &
!                          PCell%EEdge_Area(Isize,j)/PGrid%dy(Isize,j)
!      end do
      ! Compute the velocity at the boundary such that the global mass is conserved.
      do i = 1,Isize
        TVar%v(i,Jsize)=((-TVar%u(i,Jsize)*PCell%EEdge_Area(i,Jsize)           &
                +TVar%u(i-1,Jsize)*PCell%WEdge_Area(i,Jsize))*PGrid%dy(i,Jsize)&
                +TVar%v(i,Jsize-1)*PCell%SEdge_Area(i,Jsize)*PGrid%dx(i,Jsize))&
                                  /PGrid%dx(i,Jsize)/PCell%NEdge_Area(i,Jsize)
      end do
      ! Compute the mass error.   
      do i = ibeg,ibeg+Isize-1
        do j = jbeg,jbeg+Jsize-1
          TVar%mres(i,j)=(TVar%u(i,j)*PCell%EEdge_Area(i,j)-TVar%u(i-1,j)*     &
                          PCell%WEdge_Area(i,j))*PGrid%dy(i,j)                 &
                        +(TVar%v(i,j)*PCell%NEdge_Area(i,j)-TVar%v(i,j-1)*     &
                          PCell%SEdge_Area(i,j))*PGrid%dx(i,j)                 &
                         -BoomCase%vs*PCell%nyS(i,j)*PCell%WlLh(i,j)
          if(PCell%vofS(i,j)>1.d0-epsi) TVar%mres(i,j)=0.d0
        end do
      end do
      call VariablesInternalCellCondition(TVar,PCell,UCell,VCell)
      deallocate(GradPUV)
      deallocate(Pred%u)
      deallocate(Pred%v)
      deallocate(Proj%Pp)
      deallocate(PU%Dp)
      deallocate(PV%Dp)
    END SUBROUTINE UpdatePUV

    SUBROUTINE VariablesInternalCellCondition(TVar,PCell,UCell,VCell)
      IMPLICIT NONE
      TYPE(TVariables),INTENT(INOUT):: TVar
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      INTEGER(kind=it4b):: i,j
      do i = ibeg,ibeg+Isize-1
        do j = jbeg,jbeg+Jsize-1
          if(PCell%vofS(i,j)>1.d0-epsi) then
            TVar%p(i,j) = 0.d0
          end if
          if(UCell%vofS(i,j)>=1.d0-epsi) then
            TVar%u(i,j) = 0.d0
          end if
          if(VCell%vofS(i,j)>=1.d0-epsi) then
            TVar%v(i,j) = 0.d0
          end if
        end do
      end do
    END SUBROUTINE VariablesInternalCellCondition
end Module ComputePUV
