Module Mesh
    USE PrecisionVar
    USE MPI
    IMPLICIT NONE
    PRIVATE
    INTEGER(kind=it4b),PRIVATE :: ight=1,jght=1
    INTEGER(kind=it4b),PRIVATE :: ibeg=1,jbeg=1
    INTEGER(kind=it4b),PRIVATE :: Isize,Jsize

    TYPE :: TPoint
      REAL(KIND=dp)::x,y
    End TYPE TPoint

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

    Interface getMeshSizes
       Module procedure getMeshSizes
    End interface

    Interface TsimcoMesh
       Module procedure construct
    End interface

    Interface TPoint
       Module procedure construct_point
    End interface
    PUBLIC :: getMeshSizes, TPoint
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

        this%UCell=constructCell(ibeg-1,jbeg-1)
        this%VCell=constructCell(ibeg-1,jbeg-1)
        this%PCell=constructCell(0,0)
      end function construct

      TYPE(Cell) function constructCell(o1, o2) RESULT( this )
        INTEGER(it4b), INTENT(in) :: o1, o2
        allocate(this%WlLh      (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%delh      (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%vof       (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%phi       (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%nx        (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%ny        (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%vofS      (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%phiS      (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%nxS       (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%nyS       (ibeg:Isize+o1,jbeg:Jsize+o2))
        allocate(this%MoExCell  (ibeg:Isize+o1,jbeg:Jsize+o2))

        allocate(this%EEdge_Area(o1:Isize+o1+1,o2:Jsize+o2+1))
        allocate(this%WEdge_Area(o1:Isize+o1+1,o2:Jsize+o2+1))
        allocate(this%NEdge_Area(o1:Isize+o1+1,o2:Jsize+o2+1))
        allocate(this%SEdge_Area(o1:Isize+o1+1,o2:Jsize+o2+1))

        allocate(this%EtaE      (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%EtaN      (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%DAlE      (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%DAlN      (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%AlE       (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%AlN       (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%SxE       (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%SyN       (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
        allocate(this%PosNu     (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))

        allocate(this%FCE       (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
        allocate(this%FCN       (ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
        
        allocate(this%Cell_Cent (Isize,Jsize,2))
        allocate(this%MsCe      (Isize,Jsize,2))

      end function constructCell
      
      Subroutine getMeshSizes(ibege, jbege, Isizee, Jsizee, ighte, jghte)
        INTEGER(it4b), OPTIONAL,INTENT(out) :: ibege, jbege, Isizee, Jsizee, ighte, jghte
        !
        IF(PRESENT(ibege))  ibege  = ibeg
        IF(PRESENT(jbege))  jbege  = jbeg
        IF(PRESENT(Isizee)) Isizee = Isize
        IF(PRESENT(Jsizee)) Jsizee = Jsize
        IF(PRESENT(ighte))  ighte  = ight
        IF(PRESENT(jghte))  jghte  = jght
      end Subroutine getMeshSizes

      Subroutine Initial_Grid2(this, Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,   &
                                                     Lref,NonUniformMesh)
        IMPLICIT NONE
        CLASS(TsimcoMesh),INTENT(INOUT):: this
        TYPE(TPoint),INTENT(IN)  :: Start_Point,End_Point,ReS,ReE
        REAL(KIND=dp),INTENT(IN):: Lref
        INTEGER(kind=it4b),INTENT(IN):: Irec,Jrec,NI,NJ
        INTEGER(kind=it4b),INTENT(IN):: NonUniformMesh
        !
        REAL(KIND=dp),DIMENSION(:),allocatable:: x,y
        INTEGER:: i,j,IJsizeS,IJsizeE
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

      subroutine HYPRE_CreateGrid2(this)
        CLASS(TsimcoMesh),INTENT(INOUT):: this
        INTEGER:: ilower(0:1),iupper(0:1), ierr
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
      TYPE(TPoint) function construct_point(x, y) RESULT(this)
        REAL(dp), INTENT(in) :: x, y
        !
        this % x = x
        this % y = y
      end function construct_point
end module Mesh
