Module Mesh
    USE PrecisionVar
    USE MPI
    IMPLICIT NONE
    PRIVATE
    INTEGER(kind=it4b),PUBLIC::ight=1,jght=1
    INTEGER(kind=it4b),PRIVATE:: ibeg=1,jbeg=1
    INTEGER(kind=it4b),PRIVATE:: Isize,Jsize
    TYPE,PUBLIC:: Point
      REAL(KIND=dp)::x,y
    End TYPE Point

    TYPE,PUBLIC:: Grid
      INTEGER*8:: Grid_Id
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: x,y,dx,dy
      REAL(KIND=dp):: Lref
    End TYPE Grid

    TYPE,PUBLIC:: Cell
                          ! cell_TYPE: 0 fluid cell and not target cell
                          !            1 boudary cell
                          !            2 solid cell
      ! cell number for Poisson solving equation only for Pressure cell
      INTEGER(kind=it4b),DIMENSION(:,:),allocatable:: PosNu,MoExCell
      INTEGER(kind=it4b):: ExtCell
      ! Cell number of Implicit solving for diffusive term for u-velocity and v-velocity
      ! Volume fraction of liquid and level set function
      ! Volume fraction of Solid and level set function vofS,phiS
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: Vof,Phi,VofS,PhiS
      ! normal vector of Interface and Solid Boundary
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: nx,ny,nxS,nyS
      ! the center of cell
      ! the area of edge in x,y direction
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: EEdge_Area,WEdge_Area
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: NEdge_Area,SEdge_Area
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: Cell_Cent,MsCe
      ! the cellface center of velocity cell
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: WlLh,delh
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: FCE,FCN
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: EtaE,EtaN,DAlE,DAlN
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: AlE,AlN,SxE,SyN
      ! the cellface center of velocity cell
    End TYPE Cell

    TYPE, PUBLIC :: TsimcoMesh
      TYPE(Grid)    :: UGrid
      TYPE(Grid)    :: VGrid
      TYPE(Grid)    :: PGrid
      TYPE(Cell)    :: UCell
      TYPE(Cell)    :: VCell
      TYPE(Cell)    :: PCell

    CONTAINS
      PROCEDURE, PASS(this), PUBLIC :: InitialUVGrid2
      PROCEDURE, PASS(this), PUBLIC :: Initial_Grid2
      PROCEDURE, PASS(this), PUBLIC :: HYPRE_CreateGrid2

    END TYPE TsimcoMesh

!    PUBLIC:: Initial_Grid,HYPRE_CreateGrid,InitialUVGrid
!    Interface Initial_Grid
!      Module procedure Initial_Grid
!    End interface
!    Interface HYPRE_CreateGrid
!      Module procedure HYPRE_CreateGrid
!    End interface
!    Interface InitialUVGrid
!      Module procedure InitialUVGrid
!    End interface
    Interface getMeshSizes
       Module procedure getMeshSizes
    End interface

    Interface TsimcoMesh
       Module procedure construct
    End interface
    PUBLIC :: getMeshSizes
    Contains

      !Assume ibeg and jbeg are constant equal 1 for now
      TYPE(TsimcoMesh) function construct(Isize_in, Jsize_in) result ( this )
        !
        INTEGER(it4b), INTENT(in) :: Isize_in
        INTEGER(it4b), INTENT(in) :: Jsize_in
        !
        Isize = Isize_in
        Jsize = Jsize_in
        allocate(this%UGrid%x(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UGrid%y(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VGrid%x(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VGrid%y(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%PGrid%x(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%PGrid%y(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UGrid%dx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UGrid%dy(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VGrid%dx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VGrid%dy(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%PGrid%dx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%PGrid%dy(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))

        allocate(this%UCell%vof(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%phi(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%nx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%ny(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%vofS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%phiS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%nxS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%nyS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%EEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
        allocate(this%UCell%WEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
        allocate(this%UCell%NEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
        allocate(this%UCell%SEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
        allocate(this%UCell%MoExCell(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%UCell%EtaE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%EtaN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%DAlE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%DAlN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%AlE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%AlN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%SxE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%SyN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%FCE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
        allocate(this%UCell%FCN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
        allocate(this%UCell%WlLh(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
        allocate(this%UCell%delh(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
        allocate(this%UCell%PosNu(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%UCell%Cell_Cent(Isize,Jsize,2))
        allocate(this%UCell%MsCe(Isize,Jsize,2))

        allocate(this%VCell%vof(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%phi(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%nx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%ny(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%vofS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%phiS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%nxS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%nyS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%EEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
        allocate(this%VCell%WEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
        allocate(this%VCell%NEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
        allocate(this%VCell%SEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
        allocate(this%VCell%MoExCell(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
        allocate(this%VCell%EtaE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%EtaN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%DAlE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%DAlN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%AlE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%AlN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%SxE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%SyN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%FCE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
        allocate(this%VCell%FCN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
        allocate(this%VCell%WlLh(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
        allocate(this%VCell%delh(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
        allocate(this%VCell%PosNu(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%VCell%Cell_Cent(Isize,Jsize,2))
        allocate(this%VCell%MsCe(Isize,Jsize,2))

        allocate(this%PCell%vof(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%phi(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%nx(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%ny(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%vofS(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%phiS(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%nxS(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%nyS(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%EEdge_Area(ibeg-1:Isize+1,jbeg-1:Jsize+1))
        allocate(this%PCell%WEdge_Area(ibeg-1:Isize+1,jbeg-1:Jsize+1))
        allocate(this%PCell%NEdge_Area(ibeg-1:Isize+1,jbeg-1:Jsize+1))
        allocate(this%PCell%SEdge_Area(ibeg-1:Isize+1,jbeg-1:Jsize+1))
        allocate(this%PCell%FCE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
        allocate(this%PCell%FCN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
        allocate(this%PCell%WlLh(ibeg:Isize,jbeg:Jsize))
        allocate(this%PCell%PosNu(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%PCell%Cell_Cent(Isize,Jsize,2))

      end function construct
      Subroutine getMeshSizes(ibege, jbege, Isizee, Jsizee)
        INTEGER(it4b), INTENT(out) :: ibege, jbege, Isizee, Jsizee
        !
        ibege  = ibeg
        jbege  = jbeg
        Isizee = Isize
        Jsizee = Jsize
      end Subroutine getMeshSizes

!      Subroutine Initial_Grid(Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,   &
!                                                     TGrid,Lref,NonUniformMesh)
!        IMPLICIT NONE
!        TYPE(Grid),INTENT(INOUT):: TGrid
!        TYPE(Point),INTENT(IN)  :: Start_Point,End_Point,ReS,ReE
!        REAL(KIND=dp),INTENT(IN):: Lref
!        INTEGER(kind=it4b),INTENT(IN):: Irec,Jrec,NI,NJ
!        INTEGER(kind=it4b),INTENT(IN):: NonUniformMesh
!        REAL(KIND=dp),DIMENSION(:),allocatable:: x,y
!        INTEGER:: i,j,IJsizeS,IJsizeE
!        REAL(KIND=dp):: beta,dx,dy,dl
!        allocate(x(NI))
!        allocate(y(NJ))
!        TGrid%Lref = Lref
!
!      ! for y-direction
!        if(NonUniformMesh==1) then
!          dy=(ReE%y-ReS%y)/dble(Jrec-1)
!          dl=ReS%y-Start_Point%y
!       !  find the number of grid points from bottom wall to grid refining region
!          IJsizeS=(NJ-Jrec)*dabs(ReS%y-Start_Point%y)/                         &
!                             (dabs(ReS%y-Start_Point%y)+dabs(End_Point%y-ReE%y))
!          call NewtonRaphson(beta,dl,dy,IJsizeS)
!          y(jbeg)=Start_Point%y
!          do j=jbeg+1,jbeg+ibeg+IJsizeS
!            y(j)=y(j-1)+dy*beta**(IJsizeS-(j-jbeg-1))
!          end do
!       !  the grid refining region
!          do j=jbeg+IJsizeS+1,jbeg+IJsizeS+Jrec-1
!            y(j)=y(j-1)+dy
!          end do
!          dl=End_Point%y-ReE%y
!          IJsizeE=(NJ-Jrec)-IJsizeS
!       !  the section from the grid refining region to open air
!          call NewtonRaphson(beta,dl,dy,IJsizeE)
!          do j=jbeg+IJsizeS+Jrec,NJ
!            y(j)=y(j-1)+dy*beta**(j-(jbeg+IJsizeS+Jrec-1))
!          end do
!        else
!          dy=(End_Point%y-Start_Point%y)/dble(NJ-1)
!          do j=1,NJ
!            y(j)=Start_Point%y+dy*dble(j-1)
!          end do
!        end if
!
!      ! For x-direction
!        if(NonUniformMesh==1) then
!          dx=(ReE%x-ReS%x)/(Irec-1)
!          dl=ReS%x-Start_Point%x
!      !  find the number of grid points from inlet to grid refining region
!          IJsizeS=(NI-Irec)*dabs(ReS%x-Start_Point%x)/                         &
!                             (dabs(ReS%x-Start_Point%x)+dabs(End_Point%x-ReE%x))
!          call NewtonRaphson(beta,dl,dx,IJsizeS)
!          x(ibeg)=Start_Point%x
!          do i=ibeg+1,ibeg+IJsizeS
!            x(i)=x(i-1)+dx*beta**(IJsizeS-(i-ibeg-1))
!          end do
!      !  the grid refining region
!          do i=ibeg+IJsizeS+1,ibeg+IJsizeS+Irec-1
!            x(i)=x(i-1)+dx
!          end do
!          dl=End_Point%x-ReE%x
!          IJsizeE=(NI-Irec)-IJsizeS
!      !  the section from the grid refining region to outlet
!          call NewtonRaphson(beta,dl,dx,IJsizeE)
!          do i=ibeg+IJsizeS+Irec,NI
!            x(i)=x(i-1)+dx*beta**(i-(ibeg+IJsizeS+Irec-1))
!          end do
!        else
!          dx=(End_Point%x-Start_Point%x)/dble(NI-1)
!          do i=1,NI
!            x(i)=Start_Point%x+dx*dble(i-1)
!          end do
!        end if
!        do i=ibeg,ibeg+NI-2
!          do j=jbeg,jbeg+NJ-2
!            TGrid%x(i,j)=0.5d0*(x(i+1)+x(i))
!            TGrid%dx(i,j)=x(i+1)-x(i)
!            TGrid%y(i,j)=0.5d0*(y(j+1)+y(j))
!            TGrid%dy(i,j)=y(j+1)-y(j)
!          end do
!        end do
!        deallocate(x,y)
!      end subroutine Initial_Grid

      Subroutine Initial_Grid2(this, Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,   &
                                                     Lref,NonUniformMesh)
        IMPLICIT NONE
        CLASS(TsimcoMesh),INTENT(INOUT):: this
        TYPE(Point),INTENT(IN)  :: Start_Point,End_Point,ReS,ReE
        REAL(KIND=dp),INTENT(IN):: Lref
        INTEGER(kind=it4b),INTENT(IN):: Irec,Jrec,NI,NJ
        INTEGER(kind=it4b),INTENT(IN):: NonUniformMesh
        !
        REAL(KIND=dp),DIMENSION(:),allocatable:: x,y
        INTEGER:: i,j,IJsizeS,IJsizeE, ibeg, jbeg
        REAL(KIND=dp):: beta,dx,dy,dl
        allocate(x(NI))
        allocate(y(NJ))
        this%PGrid%Lref = Lref

      ! for y-direction
        if(NonUniformMesh==1) then
          dy=(ReE%y-ReS%y)/dble(Jrec-1)
          dl=ReS%y-Start_Point%y
       !  find the number of grid points from bottom wall to grid refining region
          IJsizeS=(NJ-Jrec)*dabs(ReS%y-Start_Point%y)/                         &
                             (dabs(ReS%y-Start_Point%y)+dabs(End_Point%y-ReE%y))
          call NewtonRaphson(beta,dl,dy,IJsizeS)
          y(jbeg)=Start_Point%y
          do j=jbeg+1,jbeg+ibeg+IJsizeS
            y(j)=y(j-1)+dy*beta**(IJsizeS-(j-jbeg-1))
          end do
       !  the grid refining region
          do j=jbeg+IJsizeS+1,jbeg+IJsizeS+Jrec-1
            y(j)=y(j-1)+dy
          end do
          dl=End_Point%y-ReE%y
          IJsizeE=(NJ-Jrec)-IJsizeS
       !  the section from the grid refining region to open air
          call NewtonRaphson(beta,dl,dy,IJsizeE)
          do j=jbeg+IJsizeS+Jrec,NJ
            y(j)=y(j-1)+dy*beta**(j-(jbeg+IJsizeS+Jrec-1))
          end do
        else
          dy=(End_Point%y-Start_Point%y)/dble(NJ-1)
          do j=1,NJ
            y(j)=Start_Point%y+dy*dble(j-1)
          end do
        end if

      ! For x-direction
        if(NonUniformMesh==1) then
          dx=(ReE%x-ReS%x)/(Irec-1)
          dl=ReS%x-Start_Point%x
      !  find the number of grid points from inlet to grid refining region
          IJsizeS=(NI-Irec)*dabs(ReS%x-Start_Point%x)/                         &
                             (dabs(ReS%x-Start_Point%x)+dabs(End_Point%x-ReE%x))
          call NewtonRaphson(beta,dl,dx,IJsizeS)
          x(ibeg)=Start_Point%x
          do i=ibeg+1,ibeg+IJsizeS
            x(i)=x(i-1)+dx*beta**(IJsizeS-(i-ibeg-1))
          end do
      !  the grid refining region
          do i=ibeg+IJsizeS+1,ibeg+IJsizeS+Irec-1
            x(i)=x(i-1)+dx
          end do
          dl=End_Point%x-ReE%x
          IJsizeE=(NI-Irec)-IJsizeS
      !  the section from the grid refining region to outlet
          call NewtonRaphson(beta,dl,dx,IJsizeE)
          do i=ibeg+IJsizeS+Irec,NI
            x(i)=x(i-1)+dx*beta**(i-(ibeg+IJsizeS+Irec-1))
          end do
        else
          dx=(End_Point%x-Start_Point%x)/dble(NI-1)
          do i=1,NI
            x(i)=Start_Point%x+dx*dble(i-1)
          end do
        end if
        do i=ibeg,ibeg+NI-2
          do j=jbeg,jbeg+NJ-2
            this%PGrid%x(i,j)=0.5d0*(x(i+1)+x(i))
            this%PGrid%dx(i,j)=x(i+1)-x(i)
            this%PGrid%y(i,j)=0.5d0*(y(j+1)+y(j))
            this%PGrid%dy(i,j)=y(j+1)-y(j)
          end do
        end do
        deallocate(x,y)
      end subroutine Initial_Grid2

!      subroutine InitialUVGrid(PGrid,TGrid,UV,Lref)
!        IMPLICIT NONE
!        TYPE(Grid),INTENT(IN):: PGrid
!        TYPE(Grid),INTENT(INOUT):: TGrid
!        REAL(KIND=dp),INTENT(IN):: Lref
!        INTEGER(kind=it4b),INTENT(IN):: UV
!        INTEGER:: i,j
!        ! for UGrid
!        TGrid%Lref = Lref
!        if(UV==0) then
!          do i = ibeg,Isize
!            do j = jbeg,Jsize
!              TGrid%x(i,j) = PGrid%x(i,j)+0.5d0*PGrid%dx(i,j)
!              TGrid%y(i,j) = PGrid%y(i,j)
!              TGrid%dy(i,j) = PGrid%dy(i,j)
!              if(i<Isize) then
!                TGrid%dx(i,j) = PGrid%x(i+1,j)-PGrid%x(i,j)
!              else
!                TGrid%dx(i,j) = PGrid%dx(i,j)
!              end if
!            end do
!          end do
!        else
!          do i = ibeg,Isize
!            do j = jbeg,Jsize
!              TGrid%x(i,j) = Pgrid%x(i,j)
!              TGrid%y(i,j) = PGrid%y(i,j)+0.5d0*PGrid%dy(i,j)
!              Tgrid%dx(i,j) = PGrid%dx(i,j)
!              if(j<Jsize) then
!                TGrid%dy(i,j) = PGrid%y(i,j+1)-PGrid%y(i,j)
!              else
!                TGrid%dy(i,j) = PGrid%dy(i,j)
!              end if
!            end do
!          end do
!        end if
!      end subroutine

      subroutine InitialUVGrid2(this,LrefU, LrefV)
        IMPLICIT NONE
        CLASS(TsimcoMesh),INTENT(INOUT):: this
        REAL(KIND=dp),INTENT(IN):: LrefU
        REAL(KIND=dp),INTENT(IN):: LrefV
        INTEGER:: i,j
        ! for UGrid
        this%UGrid%Lref = LrefU
        this%VGrid%Lref = LrefV
        do i = ibeg,Isize
          do j = jbeg,Jsize
            this%UGrid%x(i,j)  = this%PGrid%x(i,j)+0.5d0*this%PGrid%dx(i,j)
            this%UGrid%y(i,j)  = this%PGrid%y(i,j)
            this%UGrid%dy(i,j) = this%PGrid%dy(i,j)
            if(i<Isize) then
              this%UGrid%dx(i,j) = this%PGrid%x(i+1,j)-this%PGrid%x(i,j)
            else
              this%UGrid%dx(i,j) = this%PGrid%dx(i,j)
            end if
          end do
        end do
        do i = ibeg,Isize
          do j = jbeg,Jsize
            this%VGrid%x(i,j)  = this%Pgrid%x(i,j)
            this%VGrid%y(i,j)  = this%PGrid%y(i,j)+0.5d0*this%PGrid%dy(i,j)
            this%Vgrid%dx(i,j) = this%PGrid%dx(i,j)
            if(j<Jsize) then
              this%VGrid%dy(i,j) = this%PGrid%y(i,j+1)-this%PGrid%y(i,j)
            else
              this%VGrid%dy(i,j) = this%PGrid%dy(i,j)
            end if
          end do
        end do
      end subroutine

!      subroutine HYPRE_CreateGrid(TGrid)
!        IMPLICIT NONE
!        TYPE(Grid),INTENT(INOUT):: TGrid
!        INTEGER:: ilower(0:1),iupper(0:1)
!        ilower(0) = ibeg
!        ilower(1) = jbeg
!        iupper(0) = ibeg+Isize-1
!        iupper(1) = jbeg+Jsize-1
!
!        call HYPRE_StructGridCreate(MPI_COMM_WORLD,2,TGrid%Grid_Id,ierr)
!        call HYPRE_StructGridSetExtents(TGrid%Grid_Id,ilower,iupper,ierr)
!        call HYPRE_StructGridAssemble(TGrid%Grid_Id,ierr)
!      end subroutine HYPRE_CreateGrid

      subroutine HYPRE_CreateGrid2(this)
        CLASS(TsimcoMesh),INTENT(INOUT):: this
        INTEGER:: ilower(0:1),iupper(0:1)
        ilower(0) = ibeg
        ilower(1) = jbeg
        iupper(0) = ibeg+Isize-1
        iupper(1) = jbeg+Jsize-1

        call HYPRE_StructGridCreate(MPI_COMM_WORLD,2,this%PGrid%Grid_Id,ierr)
        call HYPRE_StructGridSetExtents(this%PGrid%Grid_Id,ilower,iupper,ierr)
        call HYPRE_StructGridAssemble(this%PGrid%Grid_Id,ierr)
      end subroutine HYPRE_CreateGrid2

      subroutine NewtonRaphson(beta,dl,dx,IJsize)
        REAL(KIND=dp),INTENT(INOUT):: beta
        REAL(KIND=dp),INTENT(IN):: dl,dx
        INTEGER(kind=it4b),INTENT(IN):: IJsize
        REAL(KIND=dp):: tol,fx,dfx
        tol = 1.d0
        beta = 1.001d0
        do while(tol>1.d-14)
          fx = dx*(beta**(IJsize+1)-beta)/(beta-1.d0)-dl
          dfx = dx*(IJsize*beta**(IJsize+1)-(IJsize+1)*beta**IJsize+1)/(beta-1)**2.d0
          tol = dabs(fx/dfx)
          beta = beta-fx/dfx
        end do
      end subroutine
end module Mesh
