Module Cutcell
    USE PrecisionVar
    USE Mesh, ONLY : TsimcoMesh, Grid, Cell, TPoint, getMeshSizes
    USE StateVariables, ONLY : TVariables, getSolverVariables
    USE Constants, ONLY : epsi, epsiF
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: Grid_Preprocess,NewCellFace
    Interface Grid_Preprocess
      Module procedure Grid_Preprocess
    End interface Grid_Preprocess
    Interface NumberExternalCell
      Module procedure NumberExternalCell
    End interface
    Interface NewCellFace
      Module procedure NewCellFace
    End interface
    Contains
      Subroutine Grid_Preprocess(simcomesh, TVar,itt)
        IMPLICIT NONE
        TYPE(TsimcoMesh), INTENT(inout) :: simcomesh
        TYPE(TVariables),INTENT(IN):: TVar
        INTEGER(kind=it8b),INTENT(IN):: itt
        INTEGER(kind=it4b):: i,j
        INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
        call getMeshSizes(ibeg, jbeg, Isize, Jsize)
        call Cell_Geo_Cal(simcomesh%PGrid,simcomesh%PCell, ibeg, jbeg, Isize, Jsize)
        call Cell_Geo_Cal(simcomesh%UGrid,simcomesh%UCell, ibeg, jbeg, Isize, Jsize)
        call Cell_Geo_Cal(simcomesh%VGrid,simcomesh%VCell, ibeg, jbeg, Isize, Jsize)
        call DefineMomentumExchangeCell(simcomesh%PCell,simcomesh%UCell,simcomesh%VCell, ibeg, &
             &                          jbeg, Isize, Jsize)
        call NumberExternalCell(simcomesh%PCell,0,0, ibeg, jbeg, Isize, Jsize)
        call NumberExternalCell(simcomesh%UCell,1,0, ibeg, jbeg, Isize, Jsize)
        call NumberExternalCell(simcomesh%VCell,0,1, ibeg, jbeg, Isize, Jsize)
      End subroutine Grid_Preprocess

      Subroutine Cell_Geo_Cal(TGrid,TCell, ibeg, jbeg, Isize, Jsize)
        IMPLICIT NONE
        TYPE(Grid),INTENT(IN):: TGrid
        TYPE(Cell),INTENT(INOUT):: TCell
        INTEGER(it4b), INTENT(in) :: ibeg, jbeg, Isize, Jsize
        TYPE(TPoint):: Pt(0:1,0:1),FaCe
        INTEGER(kind=it4b):: i,j,ii,jj,ctr
        REAL(KIND=dp):: Nodels(0:1,0:1),MaxFace,AverageArea
        REAL(KIND=dp):: dh,AverageFace,Roundoff
        REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: FCW,FCS
        LOGICAL :: ICorProb
        allocate(FCW(Isize,Jsize,2))
        allocate(FCS(Isize,Jsize,2))
        Roundoff = 1.d-24
        do i = ibeg,Isize+ibeg-1
          do j = jbeg,Jsize+jbeg-1
            TCell%Cell_Cent(i,j,1)=0.d0
            TCell%Cell_Cent(i,j,2)=0.d0
            ctr = 0
            dh = 0.d0
            do ii = 0,1
              do jj = 0,1
                Pt(ii,jj)%x=-TGrid%dx(i,j)*(0.5d0-dble(ii))
                Pt(ii,jj)%y=-TGrid%dy(i,j)*(0.5d0-dble(jj))
              end do
            end do
            call Edge_Geo_Cal(Pt(0,0),Pt(0,1),TCell%nxs(i,j),TCell%nys(i,j),   &
                                    TCell%phis(i,j),TCell%WEdge_Area(i,j),FaCe)
            FCW(i,j,1)=Face%x
            FCW(i,j,2)=Face%y
            call Edge_Geo_Cal(Pt(0,1),Pt(1,1),TCell%nxs(i,j),TCell%nys(i,j),   &
                                    TCell%phis(i,j),TCell%NEdge_Area(i,j),FaCe)
            TCell%FCN(i,j,1)=FaCe%x
            TCell%FCN(i,j,2)=FaCe%y
            call Edge_Geo_Cal(Pt(1,1),Pt(1,0),TCell%nxs(i,j),TCell%nys(i,j),   &
                                    TCell%phis(i,j),TCell%EEdge_Area(i,j),FaCe)
            TCell%FCE(i,j,1)=FaCe%x
            TCell%FCE(i,j,2)=FaCe%y
            call Edge_Geo_Cal(Pt(1,0),Pt(0,0),TCell%nxs(i,j),TCell%nys(i,j),   &
                                    TCell%phis(i,j),TCell%SEdge_Area(i,j),FaCe)
            FCS(i,j,1)=Face%x
            FCS(i,j,2)=Face%y
          ! Calculate the face area of internal cell
          ! An edge sharing with internal cell will be set up at 0
          ! The other edges will be set up at actual value.
          end do
        end do
        if(allocated(TCell%WlLh)) then
          do i=1,Isize
            do j=1,Jsize
              TCell%WlLh(i,j)=dsqrt(((TCell%EEdge_Area(i,j)-                   &
                              TCell%WEdge_Area(i,j))*TGrid%dx(i,j))**2.d0+     &
                            ((TCell%NEdge_Area(i,j)-TCell%SEdge_Area(i,j))*    &
                              TGrid%dy(i,j))**2.d0)
      !        TCell%WlLh(i,j) = 0.d0
            end do
          end do
        end if

        TCell%FCE(ibeg-1,:,1) = TCell%FCE(ibeg,:,1)-TGrid%dx(ibeg,1)
        TCell%FCE(ibeg-1,:,2) = TCell%FCE(ibeg,:,2)
        TCell%FCE(ibeg+Isize,:,1) = TCell%FCE(ibeg+Isize-1,:,1)+               &
                                                       TGrid%dx(ibeg+Isize-1,1)
        TCell%FCE(ibeg+Isize,:,2) = TCell%FCE(ibeg+Isize-1,:,2)

        TCell%FCN(ibeg-1,:,1) = TCell%FCN(ibeg,:,1)-TGrid%dx(ibeg,1)
        TCell%FCN(ibeg-1,:,2) = TCell%FCN(ibeg,:,2)
        TCell%FCN(ibeg+Isize,:,1) = TCell%FCN(ibeg+Isize-1,:,1)+               &
                                                       TGrid%dx(ibeg+Isize-1,1)
        TCell%FCN(ibeg+Isize,:,2) = TCell%FCN(ibeg+Isize-1,:,2)

        TCell%FCE(:,jbeg-1,1) = TCell%FCE(:,jbeg,1)
        TCell%FCE(:,jbeg-1,2) = TCell%FCE(:,jbeg,2)-TGrid%dy(1,jbeg)
        TCell%FCE(:,jbeg+Jsize,1) = TCell%FCE(:,jbeg+Jsize-1,1)
        TCell%FCE(:,jbeg+Jsize,2) = TCell%FCE(:,jbeg+Jsize-1,2)+               &
                                                       TGrid%dy(1,jbeg+Jsize-1)

        TCell%FCN(:,jbeg-1,1) = TCell%FCE(:,jbeg,1)
        TCell%FCN(:,jbeg-1,2) = TCell%FCE(:,jbeg,2)-TGrid%dy(1,jbeg)
        TCell%FCN(:,jbeg+Jsize,1) = TCell%FCE(:,jbeg+Jsize-1,1)
        TCell%FCN(:,jbeg+Jsize,2) = TCell%FCE(:,jbeg+Jsize-1,2)+               &
                                                       TGrid%dy(1,jbeg+Jsize-1)

        do i=1,Isize
          do j=1,Jsize
            TCell%WEdge_Area(i,j)=1.d0-TCell%WEdge_Area(i,j)
            TCell%EEdge_Area(i,j)=1.d0-TCell%EEdge_Area(i,j)
            TCell%NEdge_Area(i,j)=1.d0-TCell%NEdge_Area(i,j)
            TCell%SEdge_Area(i,j)=1.d0-TCell%SEdge_Area(i,j)
          ! Calculate the Face Ceter for Eastern Edge
          ! (Convert from solid face to fluid face)
            if(dabs(TCell%FCE(i,j,2))<epsi) then
              TCell%FCE(i,j,2)=0.d0
            else
              TCell%FCE(i,j,2)=(0.5d0-TCell%EEdge_Area(i,j)/2.d0)*             &
                                dsign(1.d0,-TCell%FCE(i,j,2))*TGrid%dy(i,j)
            end if
            if(TCell%EEdge_Area(i,j)<epsiF) TCell%FCE(i,j,2)=0.d0

          ! Calculate the Face Ceter for Western Edge
          ! (Convert from solid face to fluid face)
            if(dabs(FCW(i,j,2))<epsi) then
              FCW(i,j,2)=0.d0
            else
              FCW(i,j,2)=(0.5d0-TCell%WEdge_Area(i,j)/2.d0)*                   &
                                      dsign(1.d0,-FCW(i,j,2))*TGrid%dy(i,j)
            end if
            if(TCell%WEdge_Area(i,j)<epsiF) FCW(i,j,2)=0.d0
          ! Calculate the Face Ceter for Northern Edge
          ! (Convert from solid face to fluid face)
            if(dabs(TCell%FCN(i,j,1))<epsi) then
              TCell%FCN(i,j,1)=0.d0
            else
              TCell%FCN(i,j,1)=(0.5d0-TCell%NEdge_Area(i,j)/2.d0)*             &
                                dsign(1.d0,-TCell%FCN(i,j,1))*TGrid%dx(i,j)
            end if
            if(TCell%NEdge_Area(i,j)<epsiF) TCell%FCN(i,j,1)=0.d0
          ! Calculate the Face Ceter for Southern Edge
          ! (Convert from solid face to fluid face)
            if(dabs(FCS(i,j,1))<epsi) then
              FCS(i,j,1)=0.d0
            else
              FCS(i,j,1)=(0.5d0-TCell%SEdge_Area(i,j)/2.d0)*                   &
                                      dsign(1.d0,-FCS(i,j,1))*TGrid%dx(i,j)
            end if
            if(TCell%SEdge_Area(i,j)<epsiF) FCS(i,j,1)=0.d0
          end do
        end do

        do i=2,Isize-1
          do j=2,Jsize-1
            AverageFace=2.d0*TCell%NEdge_Area(i,j)*TCell%SEdge_Area(i,j+1)/    &
                       (TCell%NEdge_Area(i,j)+TCell%SEdge_Area(i,j+1)+Roundoff)
            TCell%FCN(i,j,1)=(TCell%FCN(i,j,1)*TCell%NEdge_Area(i,j)+          &
                            FCS(i,j+1,1)*TCell%SEdge_Area(i,j+1))/             &
                       (TCell%NEdge_Area(i,j)+TCell%SEdge_Area(i,j+1)+Roundoff)
            TCell%FCN(i,j,2)=(TCell%FCN(i,j,2)*TCell%NEdge_Area(i,j)+          &
                       (FCS(i,j+1,2)+(TGrid%y(i,j+1)-TGrid%y(i,j)))*           &
                        TCell%SEdge_Area(i,j+1))/(TCell%NEdge_Area(i,j)+       &
                        TCell%SEdge_Area(i,j+1)+Roundoff)

            if(TCell%NEdge_Area(i,j)<epsiF.and.TCell%SEdge_Area(i,j+1)<epsiF) then
              TCell%FCN(i,j,2)=TGrid%dy(i,j)/2.d0
            end if
            if(isnan(TCell%FCN(i,j,1)).or.isnan(TCell%FCN(i,j,2))) then
              print*,FCS(i,j+1,1),TCell%NEdge_Area(i,j),TCell%SEdge_Area(i,j+1)
              print*,FCS(i,j+1,2)
              print*,'Cutcell_178'
            end if
            TCell%NEdge_Area(i,j)=AverageFace
            TCell%SEdge_Area(i,j+1)=AverageFace

            AverageFace=2.d0*TCell%EEdge_Area(i,j)*TCell%WEdge_Area(i+1,j)/    &
                       (TCell%EEdge_Area(i,j)+TCell%WEdge_Area(i+1,j)+Roundoff)
            TCell%FCE(i,j,1)=(TCell%FCE(i,j,1)*TCell%EEdge_Area(i,j)+          &
                       (FCW(i+1,j,1)+(TGrid%x(i+1,j)-TGrid%x(i,j)))*           &
                        TCell%WEdge_Area(i+1,j))/(TCell%EEdge_Area(i,j)+       &
                                              TCell%WEdge_Area(i+1,j)+Roundoff)
            TCell%FCE(i,j,2)=(TCell%FCE(i,j,2)*TCell%EEdge_Area(i,j)+          &
                            FCW(i+1,j,2)*TCell%WEdge_Area(i+1,j))/             &
                       (TCell%EEdge_Area(i,j)+TCell%WEdge_Area(i+1,j)+Roundoff)
            if(TCell%EEdge_Area(i,j)<epsiF.and.TCell%WEdge_Area(i+1,j)<epsiF) then
              TCell%FCE(i,j,1)=TGrid%dx(i,j)/2.d0
            end if
            if(isnan(TCell%FCE(i,j,1)).or.isnan(TCell%FCE(i,j,2))) then
              print*,FCW(i+1,j,1),TCell%EEdge_Area(i,j),TCell%WEdge_Area(i,j+1)
              print*,FCW(i+1,j,2)
              print*,'Cutcell_198'
            end if
            TCell%EEdge_Area(i,j)=AverageFace
            TCell%WEdge_Area(i+1,j)=AverageFace
          end do
        end do
!        print*,'test cutcell 205'
!        print*,TCell%vofS(299,76)
!        print*,TCell%WEDge_Area(299,76),TCell%EEDge_Area(299,76)
!        print*,TCell%SEDge_Area(299,76),TCell%NEDge_Area(299,76)
!        print*,
!        print*,TCell%vofS(299,77)
!        print*,TCell%WEDge_Area(299,77),TCell%EEDge_Area(299,77)
!        print*,TCell%SEDge_Area(299,77),TCell%NEDge_Area(299,77)
!        print*,
   !    Modify the face area for sharing cell
        do i=1,Isize
          do j=1,Jsize
            if(TCell%vofS(i,j)>1.d0-epsi) then
              TCell%vofS(i,j)=1.d0
              TCell%WEdge_Area(i,j)=0.d0
              TCell%EEdge_ARea(i,j)=0.d0
              TCell%NEdge_ARea(i,j)=0.d0
              TCell%SEdge_ARea(i,j)=0.d0
            ! for neighbor cell
              TCell%EEdge_Area(i-1,j)=0.d0
              TCell%NEdge_Area(i,j-1)=0.d0
              TCell%WEdge_Area(i+1,j)=0.d0
              TCell%SEdge_Area(i,j+1)=0.d0
            end if
            if((TCell%WEdge_Area(i,j)+TCell%EEdge_Area(i,j)<Roundoff.or.        &
               TCell%SEdge_ARea(i,j)+TCell%NEdge_Area(i,j)<Roundoff).and.       &
               TCell%VofS(i,j)<1.d0-epsi.and.(.not.allocated(TCell%MsCe))) then
              TCell%vofS(i,j)=1.d0
              TCell%WEdge_Area(i,j)=0.d0
              TCell%EEdge_ARea(i,j)=0.d0
              TCell%NEdge_ARea(i,j)=0.d0
              TCell%SEdge_ARea(i,j)=0.d0
            ! for neighbor cell
              TCell%EEdge_Area(i-1,j)=0.d0
              TCell%NEdge_Area(i,j-1)=0.d0
              TCell%WEdge_Area(i+1,j)=0.d0
              TCell%SEdge_Area(i,j+1)=0.d0
            end if
          end do
        end do
        do i = ibeg,ibeg+Isize-1
          do j = jbeg,jbeg+Jsize-1
            if(dabs(TCell%vofS(i,j)-1.d0)<Roundoff) then
              if(TCell%WEdge_Area(i,j)==0.d0.or.TCell%EEdge_Area(i,j)==0.d0) then
                TCell%WlLh(i,j) = TGrid%dy(i,j)
                if(allocated(TCell%delh))TCell%delh(i,j) = TGrid%dx(i,j)/2.d0
              end if
              if(TCell%NEdge_ARea(i,j)==0.d0.or.TCell%SEdge_ARea(i,j)==0.d0) then
                TCell%WlLh(i,j) = TGrid%dx(i,j)
                if(allocated(TCell%delh))TCell%delh(i,j) = TGrid%dy(i,j)/2.d0
              end if
            end if
          end do
        end do
        call getSolverVariables(ICorProb_=ICorProb)
        if((ICorProb.eqv..TRUE.).and.(.not.allocated(TCell%MoExCell))) then
          do i=1,Isize
            do j=1,Jsize-4
              if(TCell%VofS(i,j)>-epsi.and.TCell%vofS(i,j)<epsi) then
                if(TCell%vofS(i,j+1)>0.6d0-epsi.and.TCell%vofS(i,j+1)<0.6d0+   &
                                                                      epsi) then
                  if(TCell%vofS(i,j+2)>0.6d0-epsi.and.TCell%vofS(i,j+2)<0.6d0+ &
                                                                      epsi) then
                    TCell%NEdge_ARea(i,j)=0.40
                    TCell%SEDGe_Area(i,j+1)=0.4d0
                    TCell%FCN(i,j,1)=-0.3d0*TGrid%dx(i,j)
                    TCell%FCN(i,j,2)=0.5d0*TGrid%dy(i,j)
                    TCell%WlLh(i,j)=0.6d0*TGrid%dx(i,j)
                    TCell%nyS(i,j)=-1.d0
                    TCell%nxS(i,j)=0.d0
                  end if
                end if
                if(TCell%vofS(i,j+1)>1.d0-epsi.and.TCell%vofS(i,j+1)<1.d0+     &
                                                                      epsi) then
                  TCell%WlLh(i,j)=TGrid%dx(i,j)
                  TCell%nyS(i,j)=-1.d0
                  TCell%nxS(i,j)=0.d0
                end if
              end if
            end do
          end do
        end if
!        print*,'test cutcell 259'
!        if((.not.allocated(TCell%MoExCell))) then
!        print*,'***********************************************'
!        print*,'this is pcell'
!        end if
!        print*,TCell%vofS(299,76)
!        print*,TCell%WEDge_Area(299,76),TCell%EEDge_Area(299,76)
!        print*,TCell%SEDge_Area(299,76),TCell%NEDge_Area(299,76)
!        print*,
    !    print*,TCell%vofS(299,77)
    !    print*,TCell%WEDge_Area(299,77),TCell%EEDge_Area(299,77)
    !    print*,TCell%SEDge_Area(299,77),TCell%NEDge_Area(299,77)
    !    pause 'end test cutcell 267'
        deallocate(FCW,FCS)
      end subroutine Cell_Geo_Cal

      SUBROUTINE Edge_Geo_Cal(Pt1,Pt2,nxx,nyy,diss,Edge_Area,FaCe)
        IMPLICIT NONE
        TYPE(Tpoint),INTENT(IN)   :: Pt1,Pt2
        TYPE(Tpoint),INTENT(OUT)  :: Face
        REAL(KIND=dp),INTENT(IN) :: nxx,nyy,diss
        REAL(KIND=dp),INTENT(OUT):: Edge_Area
        REAL(KIND=dp)            :: epsil,Lvs_Pt1,Lvs_Pt2
        Lvs_Pt1 = nxx*Pt1%x+nyy*Pt1%y+diss
        Lvs_Pt2 = nxx*Pt2%x+nyy*Pt2%y+diss
        epsil = 1.d-14
        FaCe%x = 0.5d0*(pt1%x+pt2%x)
        FaCe%y = 0.5d0*(pt1%y+pt2%y)
        if(Lvs_Pt1>=epsil.and.Lvs_Pt2>=epsil) then
          Edge_Area = 0.d0
        elseif(Lvs_Pt1<epsil.and.Lvs_Pt2<epsil) then
          Edge_Area = 1.d0
        ! the Pt2 in solid field, Pt1 in liquid field
        elseif(Lvs_Pt1>=epsil.and.Lvs_Pt2<epsil) then
          if(dabs(pt1%x-pt2%x)<=1.d-10) then
            Edge_Area=dabs(Lvs_Pt2)/(dabs(Lvs_Pt1)+dabs(Lvs_Pt2))
            FaCe%y=pt2%y-0.5d0*Edge_Area*(pt2%y-pt1%y)
          else
            Edge_Area=dabs(Lvs_Pt2)/(dabs(Lvs_Pt1)+dabs(Lvs_Pt2))
            FaCe%x=pt2%x-0.5d0*Edge_Area*(pt2%x-pt1%x)
          end if
        ! the Pt1 in solid field, Pt2 in liquid field
        elseif(Lvs_Pt1<epsil.and.Lvs_Pt2>=epsil) then
          if(dabs(pt1%x-pt2%x)<=1.d-10) then
            Edge_Area=dabs(Lvs_Pt1)/(dabs(Lvs_Pt1)+dabs(Lvs_Pt2))
            FaCe%y=pt1%y-0.5d0*Edge_Area*(pt1%y-pt2%y)
          else
            Edge_Area=dabs(Lvs_Pt1)/(dabs(Lvs_Pt1)+dabs(Lvs_Pt2))
            FaCe%x=pt1%x-0.5d0*Edge_Area*(pt1%x-pt2%x)
          end if
        end if
      end subroutine Edge_Geo_Cal

  ! numbering the pressure cell for poisson equation
      Subroutine NumberExternalCell(TCell,iu,iv, ibeg, jbeg, Isize, Jsize)
        IMPLICIT NONE
        TYPE(Cell),INTENT(INOUT):: TCell
        INTEGER(kind=it4b),INTENT(IN):: iu,iv
        INTEGER(it4b), INTENT(in) :: ibeg, jbeg, Isize, Jsize
        INTEGER(kind=it4b):: i,j
        INTEGER(kind=it4b):: ctr
        ctr = 0
        Do i = ibeg,Isize-iu
          Do j = jbeg,Jsize-iv
            if(TCell%vofS(i,j)<1.d0-epsi) then
              TCell%Posnu(i,j)=ctr
              ctr = ctr+1
            else
              TCell%Posnu(i,j)=-1
            end if
          End do
        End do
        TCell%ExtCell = ctr-1
      End subroutine NumberExternalCell

      subroutine DefineMomentumExchangeCell(PCell,UCell,VCell, ibeg, jbeg, Isize, Jsize)
        IMPLICIT NONE
        TYPE(Cell),INTENT(INOUT):: PCell
        TYPE(Cell),INTENT(INOUT):: UCell,VCell
        INTEGER(it4b), INTENT(in) :: ibeg, jbeg, Isize, Jsize
        INTEGER(kind=it4b):: i,j
        do i = ibeg,ibeg+Isize-1
          do j = jbeg,jbeg+Jsize-1
            UCell%MoExCell(i,j) = 0
            VCell%MoExCell(i,j) = 0
          ! Define the U-cell applying the momentum exchange term
            if(UCell%VofS(i,j)>=epsi.and.UCell%VofS(i,j)<1.d0-epsi) then
              if(i<ibeg+Isize-1) then
                if(PCell%VofS(i,j)>=1.d0-epsi.or.PCell%VofS(i+1,j)>=1.d0-epsi  &
                                            .or.PCell%EEdge_Area(i,j)<epsiF) then
                  UCell%MoExCell(i,j)=1
                  PCell%EEdge_Area(i,j)=0.d0
                  PCell%WEdge_Area(i+1,j)=0.d0
                end if
              else
                if(PCell%VofS(i,j)>=1.d0-epsi) UCell%MoExCell(i,j) = 1
              end if
            end if
          ! Define the V-cell applying the momentum exchange term
            if(VCell%VofS(i,j)>=epsi.and.VCell%VofS(i,j)<1.d0-epsi) then
              if(j<jbeg+Jsize-1) then
                if(Pcell%VofS(i,j)>=1.d0-epsi.or.PCell%VofS(i,j+1)>=1.d0-epsi  &
                                            .or.PCell%NEdge_Area(i,j)<epsiF) then
                  VCell%MoExCell(i,j)=1
                  PCell%NEdge_Area(i,j)=0.d0
                  PCell%SEdge_Area(i,j+1)=0.d0
                end if
              else
                if(PCell%VofS(i,j)>=1.d0-epsi) VCell%MoExCell(i,j) = 1
              end if
            end if
          end do
        end do
      end subroutine DefineMomentumExchangeCell

      Subroutine NewCellFace(simcomesh)
        IMPLICIT NONE
        TYPE(TsimcoMesh), INTENT(inout) :: simcomesh
!        TYPE(Cell) :: PCell,UCell,VCell
!        TYPE(Grid) :: PGrid,UGrid,VGrid
        INTEGER(kind=it4b):: i,j,k,ii,jj,temp, ibeg, jbeg, Isize, Jsize
        REAL(KIND=dp):: MaxFace
        INTEGER(kind=it4b),PARAMETER:: nx=5,ny=5
        REAL(KIND=dp):: dxx,dyy,xx,yy,dd,Cdis,vol,delh,Sx,Sy
        REAL(KIND=dp):: xf,yf,nxf,nyf,delh1,delh2,tol,tol1
        REAL(KIND=dp),DIMENSION(:):: phi(4)
        REAL(KIND=dp),DIMENSION(:,:):: node(6,2)
        tol = 1.d-24
        tol1 = 1.d-14
        ! This is a shortcut
!        PCell = simcomesh % PCell
!        UCell = simcomesh % UCell
!        VCell = simcomesh % VCell
!        Pgrid = simcomesh % Pgrid
!        Ugrid = simcomesh % Ugrid
!        Vgrid = simcomesh % Vgrid
        call getMeshSizes(ibeg, jbeg, Isize, Jsize)

        !
     !   dxx = simcomesh%Pgrid%delx/dble(nx)
     !   dyy = simcomesh%Pgrid%dely/dble(ny)
        temp = 0
        do i = ibeg,ibeg+Isize-1
          do j = jbeg,jbeg+Jsize-1
            if(simcomesh%Ucell%MoExCell(i,j)/=1.and.simcomesh%Ucell%VofS(i,j)>epsi) then
              simcomesh%Ucell%Cell_Cent(i,j,1)=simcomesh%Pcell%FCE(i,j,1)-simcomesh%Pgrid%dx(i,j)/2.d0
              simcomesh%Ucell%Cell_Cent(i,j,2)=simcomesh%Pcell%FCE(i,j,2)
            end if

            if(simcomesh%Vcell%MoExCell(i,j)/=1.and.simcomesh%Vcell%VofS(i,j)>epsi) then
              simcomesh%Vcell%Cell_Cent(i,j,1)=simcomesh%Pcell%FCN(i,j,1)
              simcomesh%Vcell%Cell_Cent(i,j,2)=simcomesh%Pcell%FCN(i,j,2)-simcomesh%Pgrid%dy(i,j)/2.d0
            end if
          end do
        end do
        do i=1,Isize
          do j=1,Jsize
            if(simcomesh%Ucell%MoExCell(i,j)==1) then ! for very small simcomesh%Ucell which connects to only one simcomesh%Pcell
              if(simcomesh%Ucell%vofS(i,j)<1.d0-epsi) then
                if(simcomesh%Ucell%NEdge_Area(i,j)*simcomesh%Ucell%NEdge_Area(i,j-1)==0.d0) then
                  simcomesh%Ucell%Cell_Cent(i,j,1)=simcomesh%Ucell%FCN(i,j,1)*simcomesh%Ucell%NEdge_Area(i,j)&
                                         /(simcomesh%Ucell%NEdge_Area(i,j)+tol)+      &
                                     simcomesh%Ucell%FCN(i,j-1,1)*simcomesh%Ucell%NEdge_Area(i,j-1)&
                                         /(simcomesh%Ucell%NEdge_Area(i,j-1)+tol)
                else
                  simcomesh%Ucell%Cell_Cent(i,j,2)=simcomesh%Ucell%FCE(i,j,2)*simcomesh%Ucell%EEdge_Area(i,j)&
                                         /(simcomesh%Ucell%EEdge_Area(i,j)+tol)+      &
                                     simcomesh%Ucell%FCE(i-1,j,2)*simcomesh%Ucell%EEdge_Area(i-1,j)&
                                         /(simcomesh%Ucell%EEdge_Area(i-1,j)+tol)
                  simcomesh%Ucell%Cell_Cent(i,j,1)=-(simcomesh%Ucell%nyS(i,j)*simcomesh%Ucell%Cell_Cent(i,j,2)&
                                           +simcomesh%Ucell%phiS(i,j))/(simcomesh%Ucell%nxS(i,j)+tol)
                end if
                if(simcomesh%Ucell%EEdge_Area(i,j)*simcomesh%Ucell%EEdge_Area(i-1,j)==0.d0) then
                  simcomesh%Ucell%Cell_Cent(i,j,2)=simcomesh%Ucell%FCE(i,j,2)*simcomesh%Ucell%EEdge_Area(i,j)&
                                         /(simcomesh%Ucell%EEdge_Area(i,j)+tol)+      &
                                     simcomesh%Ucell%FCE(i-1,j,2)*simcomesh%Ucell%EEdge_Area(i-1,j)&
                                         /(simcomesh%Ucell%EEdge_Area(i-1,j)+tol)
                else
                  simcomesh%Ucell%Cell_Cent(i,j,2)=-(simcomesh%Ucell%nxS(i,j)*simcomesh%Ucell%Cell_Cent(i,j,1)&
                                           +simcomesh%Ucell%phiS(i,j))/(simcomesh%Ucell%nyS(i,j)+tol)
                end if
              else
             !   simcomesh%Ucell%Cell_Cent(i,j,1) = 0.d0
             !   simcomesh%Ucell%Cell_Cent(i,j,2) = 0.d0
              end if
              call CellLinking(simcomesh%Ugrid,simcomesh%Ucell,i,j)
            end if

            if(simcomesh%Vcell%MoExCell(i,j)==1) then ! for very small simcomesh%Vcell which connects to only one simcomesh%Pcell
              if(simcomesh%Vcell%vofS(i,j)<1.d0-epsi) then
                if(simcomesh%Vcell%NEdge_Area(i,j)*simcomesh%Vcell%NEdge_Area(i,j-1)==0.d0) then
                  simcomesh%Vcell%Cell_Cent(i,j,1)=simcomesh%Vcell%FCN(i,j,1)*simcomesh%Vcell%NEdge_Area(i,j)&
                                         /(simcomesh%Vcell%NEdge_Area(i,j)+tol)+         &
                                     simcomesh%Vcell%FCN(i,j-1,1)*simcomesh%Vcell%NEdge_Area(i,j-1)&
                                         /(simcomesh%Vcell%NEdge_Area(i,j-1)+tol)
                else
                  simcomesh%Vcell%Cell_Cent(i,j,2)=simcomesh%Vcell%FCE(i,j,2)*simcomesh%Vcell%EEdge_Area(i,j)&
                                         /(simcomesh%Vcell%EEdge_Area(i,j)+tol)+         &
                                     simcomesh%Vcell%FCE(i-1,j,2)*simcomesh%Vcell%EEdge_Area(i-1,j)&
                                         /(simcomesh%Vcell%EEdge_Area(i-1,j)+tol)
                  simcomesh%Vcell%Cell_Cent(i,j,1)=-(simcomesh%Vcell%Cell_Cent(i,j,2)*simcomesh%Vcell%nyS(i,j)&
                                           +simcomesh%Vcell%phiS(i,j))/(simcomesh%Vcell%nxS(i,j)+tol)

                end if
                if(simcomesh%Vcell%EEdge_Area(i,j)*simcomesh%Vcell%EEdge_Area(i-1,j)==0.d0) then
                  simcomesh%Vcell%Cell_Cent(i,j,2)=simcomesh%Vcell%FCE(i,j,2)*simcomesh%Vcell%EEdge_Area(i,j)&
                                         /(simcomesh%Vcell%EEdge_Area(i,j)+tol)+         &
                                     simcomesh%Vcell%FCE(i-1,j,2)*simcomesh%Vcell%EEdge_Area(i-1,j)&
                                         /(simcomesh%Vcell%EEdge_Area(i-1,j)+tol)
                else
                  simcomesh%Vcell%Cell_Cent(i,j,2)=-(simcomesh%Vcell%Cell_Cent(i,j,1)*simcomesh%Vcell%nxS(i,j)&
                                           +simcomesh%Vcell%phiS(i,j))/(simcomesh%Vcell%nyS(i,j)+tol)
                end if
              else
            !    simcomesh%Vcell%Cell_Cent(i,j,1) = 0.d0
            !    simcomesh%Vcell%Cell_Cent(i,j,2) = 0.d0
              end if
              call CellLinking(simcomesh%Vgrid,simcomesh%Vcell,i,j)
            end if
          end do
        end do
        do i = 1,Isize
          do j = 1,Jsize
          ! For simcomesh%Ucell
        !    if(i<Isize) then
        !      delh1=(simcomesh%Ucell%Cell_Cent(i,j,1)+simcomesh%Pgrid%dx(i,j)/2.d0)*simcomesh%Pcell%nxS(i,j)+&
        !             simcomesh%Ucell%Cell_Cent(i,j,2)*simcomesh%Pcell%nyS(i,j)+simcomesh%Pcell%phiS(i,j)
        !     delh2=(simcomesh%Ucell%Cell_Cent(i,j,1)-simcomesh%Pgrid%dx(i+1,j)/2.d0)*             &
        !             simcomesh%Pcell%nxS(i+1,j)+simcomesh%Ucell%Cell_Cent(i,j,2)*simcomesh%Pcell%nyS(i+1,j)+ &
        !                                                      simcomesh%Pcell%phiS(i+1,j)
        !      simcomesh%Ucell%delh(i,j)=0.5d0*dabs(delh1+delh2)+tol1
        !      if(simcomesh%Pcell%vofS(i,j)>1.d0-epsi.or.simcomesh%Pcell%vofS(i,j)<epsi) then
        !        simcomesh%Ucell%delh(i,j)=dabs(delh2)+tol1
        !      end if
        !      if(simcomesh%Pcell%vofS(i+1,j)>1.d0-epsi.or.simcomesh%Pcell%vofS(i+1,j)<epsi) then
        !        simcomesh%Ucell%delh(i,j)=dabs(delh1)+tol1
        !      end if
        !    else
              simcomesh%Ucell%delh(i,j)=dabs(simcomesh%Ucell%Cell_Cent(i,j,1)*simcomesh%Ucell%nxS(i,j)+      &
                    simcomesh%Ucell%Cell_Cent(i,j,2)*simcomesh%Ucell%nyS(i,j)+simcomesh%Ucell%phiS(i,j))+tol
              if(simcomesh%Ucell%MoExCell(i,j)==1) simcomesh%Ucell%delh(i,j)=1.d-20
        !    end if
          ! For simcomesh%Vcell
        !    if(j<Jsize) then
        !      delh1=simcomesh%Vcell%Cell_Cent(i,j,1)*simcomesh%Pcell%nxS(i,j)+                     &
        !           (simcomesh%Vcell%Cell_Cent(i,j,2)+simcomesh%Pgrid%dy(i,j)/2.d0)*                &
        !                                   simcomesh%Pcell%nyS(i,j)+simcomesh%Pcell%phiS(i,j)
        !      delh2=simcomesh%Vcell%Cell_Cent(i,j,1)*simcomesh%Pcell%nxS(i,j+1)+                   &
        !           (simcomesh%Vcell%Cell_Cent(i,j,2)-simcomesh%Pgrid%dy(i,j+1)/2.d0)*              &
        !                                   simcomesh%Pcell%nyS(i,j+1)+simcomesh%Pcell%phiS(i,j+1)
        !      simcomesh%Vcell%delh(i,j)=0.5d0*dabs(delh1+delh2)+tol1
        !      if(simcomesh%Pcell%vofS(i,j)>1.d0-epsi.or.simcomesh%Pcell%vofS(i,j)<epsi) then
       !         simcomesh%Vcell%delh(i,j)=dabs(delh2)+tol1
       !       end if
       !       if(simcomesh%Pcell%vofS(i,j+1)>1.d0-epsi.or.simcomesh%Pcell%vofS(i,j+1)<epsi) then
       !         simcomesh%Vcell%delh(i,j)=dabs(delh1)+tol1
       !       end if
       !     else
              simcomesh%Vcell%delh(i,j)=dabs(simcomesh%Vcell%Cell_Cent(i,j,1)*simcomesh%Vcell%nxS(i,j)+      &
                    simcomesh%Vcell%Cell_Cent(i,j,2)*simcomesh%Vcell%nyS(i,j)+simcomesh%Vcell%phiS(i,j))+tol
              if(simcomesh%Vcell%MoExCell(i,j)==1) simcomesh%Vcell%delh(i,j)=1.d-20
        !    end if
            simcomesh%Ucell%WlLh(i,j)=dsqrt(((simcomesh%Ucell%EEdge_Area(i,j)-simcomesh%Ucell%WEdge_Area(i,j))&
                                 *simcomesh%Ugrid%dx(i,j))**2.d0+((simcomesh%Ucell%NEdge_Area(i,j)-&
                                  simcomesh%Ucell%SEdge_Area(i,j))*simcomesh%Ugrid%dy(i,j))**2.d0)
            simcomesh%Vcell%WlLh(i,j)=dsqrt(((simcomesh%Vcell%EEdge_Area(i,j)-simcomesh%Vcell%WEdge_Area(i,j))&
                                 *simcomesh%Vgrid%dx(i,j))**2.d0+((simcomesh%Vcell%NEdge_Area(i,j)-&
                                  simcomesh%Vcell%SEdge_Area(i,j))*simcomesh%Vgrid%dy(i,j))**2.d0)
          end do
        end do
      ! Define other Coefficients which are used for convective flux and diffusive calculation
      ! For simcomesh%Ucell
        call EastFaceInterpolationInf(simcomesh%Ucell,simcomesh%Ugrid,simcomesh%Pgrid,1, Isize, Jsize)
        call NorthFaceInterpolationInf(simcomesh%Ucell,simcomesh%Ugrid,simcomesh%Vgrid,0, Isize, Jsize)
        call EastFaceInterpolationInf(simcomesh%Vcell,simcomesh%Vgrid,simcomesh%Ugrid,0, Isize, Jsize)
        call NorthFaceInterpolationInf(simcomesh%Vcell,simcomesh%Vgrid,simcomesh%Pgrid,1, Isize, Jsize)
      end subroutine NewCellFace

      SUBROUTINE EastFaceInterpolationInf(TCell,TGrid,BGrid,iu, Isize, Jsize)
        IMPLICIT NONE
        TYPE(Grid),INTENT(IN):: TGrid,BGrid
        TYPE(Cell),INTENT(INOUT):: TCell
        INTEGER(kind=it4b),INTENT(IN):: iu
        INTEGER(it4b), INTENT(in) :: Isize, Jsize
        INTEGER(kind=it4b):: i,j
        REAL(KIND=dp):: IntPointDist,FaceCenterDist,Sy,xf,yf,nxf,nyf,Sx
        do j = 1,Jsize
          do i = 1,Isize-1
            if(TCell%Posnu(i,j)/=-1.and.TCell%Posnu(i+1,j)/=-1) then
              TCell%SxE(i,j)=TCell%Cell_Cent(i+1,j,1)+BGrid%dx(i+iu,j)-        &
                                                        TCell%Cell_Cent(i,j,1)
              Sy=TCell%Cell_Cent(i+1,j,2)-TCell%Cell_Cent(i,j,2)
              TCell%EtaE(i,j)=dabs(TCell%FCE(i,j,1)-TCell%Cell_Cent(i,j,1))/   &
                                                                 TCell%SxE(i,j)
              if(TCell%EtaE(i,j)>=1.d0) TCell%EtaE(i,j)=0.5d0
              xf=(1.d0-TCell%EtaE(i,j))*TCell%Cell_Cent(i,j,1)+TCell%EtaE(i,j)*&
                                    (TCell%Cell_Cent(i+1,j,1)+BGrid%dx(i+iu,j))
              yf=(1.d0-TCell%EtaE(i,j))*TCell%Cell_Cent(i,j,2)+TCell%EtaE(i,j)*&
                                     TCell%Cell_Cent(i+1,j,2)
        !      nxf=xf/dsqrt(xf**2.d0+yf**2.d0)
        !      nyf=yf/dsqrt(xf**2.d0+yf**2.d0)
              if(dabs(Sy)>0.0001d0*TGrid%dy(i,j)) then
                nxf=0.5d0*(TCell%nxS(i,j)+TCell%nxS(i+1,j))
                nyf=0.5d0*(TCell%nyS(i,j)+TCell%nyS(i+1,j))
                IntPointDist=0.5d0*(dabs(xf*TCell%nxS(i,j)+yf*TCell%nyS(i,j)+  &
                    TCell%phiS(i,j))+dabs((xf-BGrid%dx(i+iu,j))*               &
                    TCell%nxS(i+1,j)+yf*TCell%nyS(i+1,j)+TCell%phiS(i+1,j)))
                FaceCenterDist=0.5d0*(dabs(TCell%FCE(i,j,1)*TCell%nxS(i,j)+    &
                    TCell%FCE(i,j,2)*TCell%nyS(i,j)+TCell%phiS(i,j))+          &
                    dabs((TCell%FCE(i,j,1)-BGrid%dx(i+iu,j))*TCell%nxS(i+1,j)+ &
                    TCell%FCE(i,j,2)*TCell%nyS(i+1,j)+TCell%phiS(i+1,j)))
                TCell%AlE(i,j)=dabs(FaceCenterDist)/dabs(IntPointDist)
                TCell%DAlE(i,j)=Sy*nyf/(TCell%SxE(i,j)*dabs(IntPointDist))
  !              if(TCell%AlE(i,j)>1.5d0) TCell%AlE(i,j)=1.5d0
                if(dabs(IntPointDist)<1.d-5) then
                  TCell%DalE(i,j)=0.d0
                  TCell%AlE(i,j)=1.d0
                end if
                if(TCell%AlE(i,j)>=2.d0) TCell%AlE(i,j)=2.d0
              !  if(TCell%MoExCell(i,j)==1.or.TCell%MoExCell(i+1,j)==1) then
              !    TCell%DAlE(i,j) = 0.d0
              !    TCell%AlE(i,j) = 1.d0
              !  end if
              else
                TCell%DAlE(i,j)=0.d0
                TCell%AlE(i,j)=1.d0
              end if
         !     TCell%AlE(i,j)=1.d0
            else
              TCell%AlE(i,j)=1.d0
              TCell%DalE(i,j)=0.d0
              TCell%SxE(i,j)=BGrid%dx(i+iu,j)
              TCell%EtaE(i,j)=0.5d0
            end if
          end do
          TCell%SxE(i,j) = TCell%SxE(i-1,j)
          TCell%EtaE(i,j) = TCell%EtaE(i-1,j)
          TCell%DAlE(i,j) = TCell%DAlE(i-1,j)
          TCell%AlE(i,j) = TCell%AlE(i-1,j)
          TCell%SxE(0,j) = TCell%SxE(1,j)
          TCell%EtaE(0,j) = TCell%EtaE(1,j)
          TCell%DAlE(0,j) = TCell%DAlE(1,j)
          TCell%AlE(0,j) = TCell%AlN(1,j)
        end do
      END SUBROUTINE EastFaceInterpolationInf

      SUBROUTINE NorthFaceInterpolationInf(TCell,TGrid,BGrid,iv, Isize, Jsize)
        IMPLICIT NONE
        TYPE(Grid),INTENT(IN):: TGrid,BGrid
        TYPE(Cell),INTENT(INOUT):: TCell
        INTEGER(kind=it4b),INTENT(IN):: iv
        INTEGER(it4b), INTENT(in) :: Isize, Jsize
        INTEGER(kind=it4b):: i,j
        REAL(KIND=dp):: IntPointDist,FaceCenterDist,xf,yf,nxf,nyf,Sx
        do i = 1,Isize
          do j = 1,Jsize-1
            if(TCell%Posnu(i,j)/=-1.and.TCell%Posnu(i,j+1)/=-1) then
              TCell%SyN(i,j)=TCell%Cell_Cent(i,j+1,2)+BGrid%dy(i,j+iv)-        &
                                                      TCell%Cell_Cent(i,j,2)
              Sx=TCell%Cell_Cent(i,j+1,1)-TCell%Cell_Cent(i,j,1)
              TCell%EtaN(i,j)=dabs(TCell%FCN(i,j,2)-TCell%Cell_Cent(i,j,2))/   &
                                                    TCell%SyN(i,j)
              if(TCell%EtaN(i,j)>=1.d0) TCell%EtaN(i,j) = 0.5d0
              xf=(1.d0-TCell%EtaN(i,j))*TCell%Cell_Cent(i,j,1)+TCell%EtaN(i,j)*&
                                    TCell%Cell_Cent(i,j+1,1)
              yf=(1.d0-TCell%EtaN(i,j))*TCell%Cell_Cent(i,j,2)+TCell%EtaN(i,j)*&
                                    (TCell%Cell_Cent(i,j+1,2)+BGrid%dy(i,j+iv))
          !    nxf=xf/dsqrt(xf**2.d0+yf**2.d0)
          !    nyf=yf/dsqrt(xf**2.d0+yf**2.d0)
              if(dabs(Sx)>=0.0001d0*TGrid%dx(i,j)) then
                nxf=0.5d0*(TCell%nxS(i,j)+TCell%nxS(i,j+1))
                nyf=0.5d0*(TCell%nyS(i,j)+TCell%nyS(i,j+1))
                IntPointDist=0.5d0*(dabs(xf*TCell%nxS(i,j)+yf*TCell%nyS(i,j)+  &
                   TCell%phiS(i,j))+dabs(xf*TCell%nxS(i,j+1)+                  &
                   (yf-BGrid%dy(i,j+iv))*TCell%nyS(i,j+1)+TCell%phiS(i,j+1)))
                FaceCenterDist=0.5d0*(dabs(TCell%FCN(i,j,1)*TCell%nxS(i,j)+    &
                   TCell%FCN(i,j,2)*TCell%nyS(i,j)+TCell%phiS(i,j))+           &
                   dabs(TCell%FCN(i,j,1)*TCell%nxS(i,j+1)+(TCell%FCN(i,j,2)-   &
                   BGrid%dy(i,j+iv))*TCell%nyS(i,j+1)+TCell%phiS(i,j+1)))
                TCell%AlN(i,j)=dabs(IntPointDist)/dabs(FaceCenterDist)
                TCell%DAlN(i,j)=Sx*nxf/(TCell%SyN(i,j)*dabs(IntPointDist))
           !    if(TCell%AlN(i,j)>1.5d0) TCell%AlN=1.5d0
                if(dabs(IntPointDist)<1.d-5) then
                  TCell%DAlN(i,j)=0.d0
                  TCell%AlN(i,j)=1.d0
                end if
                if(TCell%AlN(i,j)>=2.d0) TCell%AlN(i,j)=2.d0
               ! if(TCell%MoExCell(i,j)==1.or.TCell%MoExCell(i,j+1)==1) then
               !   TCell%DAlN(i,j)=0.d0
               !   TCell%AlN(i,j)=1.d0
               ! end if
              else
                TCell%DAlN(i,j)=0.d0
                TCell%AlN(i,j)=1.d0
              end if
              if(isnan(TCell%DAlN(i,j))) then
                print*,i,j
                print*,iv
                print*, ""
                print*,TCell%MoExCell(i,j)
                print*, ""
                print*,TCell%nxS(i,j),TCell%nyS(i,j)
                print*,TCell%phiS(i,j)
                print*, ""
                print*,TCell%vofS(i,j),TCell%vofS(i,j+1)
                print*, ""
                print*,TCell%Cell_Cent(i,j+1,1),TCell%Cell_Cent(i,j,1)
                print*,TCell%Cell_Cent(i,j+1,2),TCell%Cell_Cent(i,j,2)
                print*, ""
                print*,Sx,nxf
                print*,'CutCell_615'
              end if
        !      TCell%AlN(i,j)=1.d0
            else
              TCell%EtaN(i,j)=0.5d0
              TCell%AlN(i,j)=1.d0
              TCell%DAlN(i,j)=0.d0
              TCell%SyN(i,j)=BGrid%dy(i,j+iv)
            end if
          end do
          TCell%SyN(i,j)=TCell%SyN(i,j-1)
          TCell%EtaN(i,j)=TCell%EtaN(i,j-1)
          TCell%DAlN(i,j)=TCell%DAlN(i,j-1)
          TCell%AlN(i,j)=TCell%AlN(i,j-1)
          TCell%SyN(i,0)=TCell%SyN(i,1)
          TCell%EtaN(i,0)=TCell%EtaN(i,1)
          TCell%DAlN(i,0)=TCell%DAlN(i,1)
          TCell%AlN(i,0)=TCell%AlN(i,1)
        end do
      END SUBROUTINE NorthFaceInterpolationInf

      SUBROUTINE CellLinking(TGrid,TCell,i,j)
        IMPLICIT NONE
        INTEGER(kind=it4b),INTENT(IN):: i,j
        TYPE(Grid),INTENT(IN):: TGrid
        TYPE(Cell),INTENT(INOUT):: TCell
        INTEGER(kind=it4b):: ii,jj
        REAL(KIND=dp):: tol
        tol=1.d-10
        ii=i+INT(sign(1.d0,TCell%nxS(i,j)))
        jj=j+INT(sign(1.d0,TCell%nyS(i,j)))
        if(TCell%vofS(ii,j)<TCell%VofS(i,jj)) then
          TCell%MsCe(i,j,1)=ii
          TCell%MsCe(i,j,2)=j
       !   TCell%Cell_Cent(i,j,1)=TCell%Cell_Cent(ii,j,1)+                      &
       !    dabs(TGrid%x(ii,j)-TGrid%x(i,j))*(1.d0-tol)*dsign(1.d0,TCell%nxS(i,j))
       !   TCell%Cell_Cent(i,j,2)=TCell%Cell_Cent(ii,j,2)
        else
          TCell%MsCe(i,j,1)=i
          TCell%MsCe(i,j,2)=jj
       !   TCell%Cell_Cent(i,j,1)=TCell%Cell_Cent(i,jj,1)
       !   TCell%Cell_Cent(i,j,2)=TCell%Cell_Cent(i,jj,2)+                      &
       !    dabs(TGrid%y(i,jj)-TGrid%y(i,j))*(1.d0-tol)*dsign(1.d0,TCell%nyS(i,j))
        end if
      END SUBROUTINE CellLinking
End module Cutcell

