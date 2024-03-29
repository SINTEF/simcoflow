Module Mesh
    USE PrecisionVar
    USE MPI
    PRIVATE
    INTEGER(kind=it4b),PUBLIC:: ibeg=1,jbeg=1
    INTEGER(kind=it4b),PUBLIC:: Isize,Jsize
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
    End TYPE

    PUBLIC:: Initial_Grid,HYPRE_CreateGrid,InitialUVGrid
    Interface Initial_Grid
      Module procedure Initial_Grid
    End interface
    Interface HYPRE_CreateGrid
      Module procedure HYPRE_CreateGrid
    End interface
    Interface InitialUVGrid
      Module procedure InitialUVGrid
    End interface
    Contains
      Subroutine Initial_Grid(Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,   &
                                                     TGrid,Lref,NonUniformMesh)
        IMPLICIT NONE
        TYPE(Grid),INTENT(INOUT):: TGrid
        TYPE(Point),INTENT(IN)  :: Start_Point,End_Point,ReS,ReE
        REAL(KIND=dp),INTENT(IN):: Lref
        INTEGER(kind=it4b),INTENT(IN):: Irec,Jrec,NI,NJ
        INTEGER(kind=it4b),INTENT(IN):: NonUniformMesh
        REAL(KIND=dp),DIMENSION(:),allocatable:: x,y
        INTEGER:: i,j,IJsizeS,IJsizeE
        REAL(KIND=dp):: beta,dx,dy,dl
        allocate(x(NI))
        allocate(y(NJ))
        TGrid%Lref = Lref

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
            TGrid%x(i,j)=0.5d0*(x(i+1)+x(i))
            TGrid%dx(i,j)=x(i+1)-x(i)
            TGrid%y(i,j)=0.5d0*(y(j+1)+y(j))
            TGrid%dy(i,j)=y(j+1)-y(j)
          end do
        end do
        deallocate(x,y)
      end subroutine Initial_Grid

      subroutine InitialUVGrid(PGrid,TGrid,UV,Lref)
        IMPLICIT NONE
        TYPE(Grid),INTENT(IN):: PGrid
        TYPE(Grid),INTENT(INOUT):: TGrid
        REAL(KIND=dp),INTENT(IN):: Lref
        INTEGER(kind=it4b),INTENT(IN):: UV
        INTEGER:: i,j
        ! for UGrid
        TGrid%Lref = Lref
        if(UV==0) then
          do i = ibeg,Isize
            do j = jbeg,Jsize
              TGrid%x(i,j) = PGrid%x(i,j)+0.5d0*PGrid%dx(i,j)
              TGrid%y(i,j) = PGrid%y(i,j)
              TGrid%dy(i,j) = PGrid%dy(i,j)
              if(i<Isize) then
                TGrid%dx(i,j) = PGrid%x(i+1,j)-PGrid%x(i,j)
              else
                TGrid%dx(i,j) = PGrid%dx(i,j)
              end if
            end do
          end do
        else
          do i = ibeg,Isize
            do j = jbeg,Jsize
              TGrid%x(i,j) = Pgrid%x(i,j)
              TGrid%y(i,j) = PGrid%y(i,j)+0.5d0*PGrid%dy(i,j)
              Tgrid%dx(i,j) = PGrid%dx(i,j)
              if(j<Jsize) then
                TGrid%dy(i,j) = PGrid%y(i,j+1)-PGrid%y(i,j)
              else
                TGrid%dy(i,j) = PGrid%dy(i,j)
              end if
            end do
          end do
        end if
      end subroutine

      subroutine HYPRE_CreateGrid(TGrid)
        IMPLICIT NONE
        TYPE(Grid),INTENT(INOUT):: TGrid
        INTEGER:: ilower(0:1),iupper(0:1)
        ilower(0) = ibeg
        ilower(1) = jbeg
        iupper(0) = ibeg+Isize-1
        iupper(1) = jbeg+Jsize-1

        call HYPRE_StructGridCreate(MPI_COMM_WORLD,2,TGrid%Grid_Id,ierr)
        call HYPRE_StructGridSetExtents(TGrid%Grid_Id,ilower,iupper,ierr)
        call HYPRE_StructGridAssemble(TGrid%Grid_Id,ierr)
      end subroutine HYPRE_CreateGrid

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
