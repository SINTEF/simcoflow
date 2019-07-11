Module PrintResult
    USE PrecisionVar
    USE Mesh, ONLY : Grid, Cell, getMeshSizes
    USE StateVariables, ONLY : TVariables
    USE Constants, ONLY : pi, Ktw, roa, row
    USE VTR, ONLY : VTR_file_handle, VTR_open_file, VTR_write_var, VTR_write_mesh, VTR_close_file
    USE Particles, ONLY : TParticle
    IMPLICIT NONE
    PRIVATE
    CHARACTER*70,PRIVATE           ::dir
    character(len=1), PARAMETER :: newline=achar(10)
    PUBLIC:: Print_Result_Tecplot_PCent,Print_Result_Tecplot_UCent,            &
             Print_Result_Tecplot_VCent,Print_Result_VTK_2D,ReadOldDataPCell,  &
             ReadOldDataVelocityCell,ReadOldDataParticle, setDir, getDir
    Interface Print_Result_Tecplot_PCent
       Module procedure Print_Result_Tecplot_PCent
    End Interface Print_Result_Tecplot_PCent
    Interface Print_Result_Tecplot_UCent
       Module procedure Print_Result_Tecplot_UCent
    End Interface Print_Result_Tecplot_UCent
    Interface Print_Result_Tecplot_VCent
       Module procedure Print_Result_Tecplot_VCent
    End Interface Print_Result_Tecplot_VCent
    Interface Print_Result_VTK_2D
       Module procedure Print_Result_VTK_2D
    End interface
    Interface ReadOldDataPCell
       Module procedure ReadOldDataPCell
    End interface
    Interface ReadOldDataVelocityCell
       Module procedure ReadOldDataVelocityCell
    End interface
    Interface ReadOldDataParticle
       Module procedure ReadOldDataParticle
    End interface
    Interface getDir
       Module procedure getDir
    End interface getDir
    Interface setDir
       Module procedure setDir
    End interface setDir
    Contains
    Subroutine Print_Result_Tecplot_PCent(TGrid,TVar,TCell,TraPar,FluxP,iter,PriPar)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: TGrid
      TYPE(TVariables),INTENT(INOUT):: TVar
      TYPE(Cell),INTENT(IN):: TCell
      TYPE(TParticle),INTENT(IN):: TraPar
      REAL(KIND=dp),INTENT(IN),DIMENSION(:,:,:),ALLOCATABLE:: FluxP
      INTEGER(kind=it8b),INTENT(IN):: iter
      INTEGER(kind=it4b),INTENT(IN):: PriPar
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: p
      Character(15):: curd
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      Allocate(p(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
      Do i = ibeg,ibeg+Isize-1
        Do j = jbeg,jbeg+Jsize-1
          If(TCell%vof(i,j)>=0.5d0) then
            p(i,j) = TVar%p(i,j)!+row/Roref*1.d0/Fr**2.d0*(End_Point%y-TGrid%y(i,j))
          Else
            p(i,j) = TVar%p(i,j)!+roa/Roref*1.d0/Fr**2.d0*(End_Point%y-TGrid%y(i,j))
          End if
        End do
      End do
      Write(curd,'(i8.8)') iter
      Open(unit=5,file=trim(adjustl(dir))//'Tecplot/Pressure_'//              &
                                             trim(curd)//'.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","p","u","v","phi","vof","phiS","vofS"',   &
                 ',"vofA","Mres","FluxU","FluxV","nxS","nyS"'
      Write(5,1112) Isize,Jsize
      Write(5,"(f24.14)")((TGrid%x(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TGrid%y(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((p(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      write(5,"(f24.14)")((0.5d0*(TVar%u(i-1,j)+TVar%u(i,j)),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      write(5,"(f24.14)")((0.5d0*(TVar%v(i,j-1)+TVar%v(i,j)),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%phi(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%vof(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%phiS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%vofS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((1.d0-TCell%vof(i,j)-TCell%vofS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TVar%mres(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      write(5,"(f24.14)")((FluxP(i,j,1),i=ibeg,Isize),j=jbeg,Jsize)
      write(5,"(f24.14)")((FluxP(i,j,2),i=ibeg,Isize),j=jbeg,Jsize)
      write(5,"(f24.14)")((TCell%nxS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      write(5,"(f24.14)")((TCell%nyS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Close(5)
1112  format(1x,'zone i=',i5,',k=',i5,',f=block')
   !  For Paraview
      if(PriPar==1) then
        Open(unit=5,file=trim(adjustl(dir))//'Paraview/ParticlesVTR_'//       &
                                             trim(curd)//'.txt',action='write')
        Write(5,*)'x,y,dp,up,vp'
        do i=1,Trapar%np
          write(5,70) TraPar%Posp(i)%x/Tgrid%Lref,',',TraPar%Posp(i)%y/Tgrid%Lref,',',             &
          TraPar%dp(i)/Tgrid%Lref,',',TraPar%uvp(i)%u/TVar%Uref,',',TraPar%uvp(i)%v/TVar%Uref
        enddo
        close(5)
      end if
   !  For Tecplot
      if(PriPar==1) then
        Open(unit=5,file=trim(adjustl(dir))//'Tecplot/'//                       &
                                                'Particles_'//trim(curd)//'.dat',action='write')
        Write(5,*)'variables = "x","y","dp","up","vp"'
        Write(5,*)'zone f=block, i=',TraPar%np,'j=',1
        Write(5,"(f24.14)")(TraPar%Posp(i)%x/Tgrid%Lref,i=1,TraPar%np)
        Write(5,"(f24.14)")(TraPar%Posp(i)%y/Tgrid%Lref,i=1,TraPar%np)
        Write(5,"(f24.14)")(TraPar%dp(i)/Tgrid%Lref,i=1,TraPar%np)
        Write(5,"(f24.14)")(TraPar%uvp(i)%u/TVar%Uref,i=1,TraPar%np)
        Write(5,"(f24.14)")(TraPar%uvp(i)%v/TVar%Uref,i=1,TraPar%np)
        Close(5)
      end if
      Deallocate(p)
 70   format((f24.14),(a3),(f24.14),(a3),(f24.14),(a3),(f24.14),(a3),(f24.14))
    End Subroutine Print_Result_Tecplot_PCent

    Subroutine Print_Result_Tecplot_UCent(TGrid,TVar,TCell,itt)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: TGrid
      TYPE(TVariables),INTENT(IN):: TVar
      TYPE(Cell),INTENT(IN):: TCell
      INTEGER(kind=it8b),INTENT(IN):: itt
      INTEGER(kind=it4b) i,j
      Character(15) curd
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      Write(curd,'(i8.8)') itt
      Open(unit=5,file=trim(adjustl(dir))//'Tecplot/'//                         &
                    'Uvelocity_'//trim(curd)//'.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","u","phi","vof","phiS","vofS","vofA"'
      Write(5,1113) Isize,Jsize
      Write(5,"(f24.14)")((TGrid%x(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TGrid%y(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TVar%u(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%phi(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%vof(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%phiS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%vofS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((1.d0-TCell%vof(i,j)-TCell%vofS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Close(5)
1113 format(1x,'zone i=',i5,',k=',i5,',f=block')
    End Subroutine Print_Result_Tecplot_UCent

    Subroutine Print_Result_Tecplot_VCent(TGrid,TVar,TCell,itt)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: TGrid
      TYPE(TVariables),INTENT(IN):: TVar
      TYPE(Cell),INTENT(IN):: TCell
      INTEGER(kind=it8b),INTENT(IN):: itt
      INTEGER(kind=it4b) i,j
      Character(15) curd
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      Write(curd,'(i8.8)') itt
      Open(unit=5,file=trim(adjustl(dir))//'Tecplot/'//                        &
                  'Vvelocity_'//trim(curd)//'.dat',action='write')
      Write(5,*) 'title = "flow"'
      Write(5,*) 'variables ="x","y","v","phi","vof","phiS","vofS","vofA"'
      Write(5,1114) Isize,Jsize
      Write(5,"(f24.14)")((TGrid%x(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TGrid%y(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TVar%v(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%phi(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%vof(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%phiS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((TCell%vofS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Write(5,"(f24.14)")((1.d0-TCell%vof(i,j)-TCell%vofS(i,j),i=ibeg,Isize+ibeg-1),j=jbeg,Jsize+jbeg-1)
      Close(5)
1114 format(1x,'zone i=',i5,',k=',i5,',f=block')
    End Subroutine Print_Result_Tecplot_VCent

    subroutine Print_Result_VTK_2D(TGrid,TVar,TCell,itt)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: TGrid
      TYPE(TVariables),INTENT(IN):: TVar
      TYPE(Cell),INTENT(IN):: TCell
      INTEGER(kind=it8b),INTENT(IN):: itt
      TYPE(VTR_file_handle):: fd
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      call VTR_open_file(Prefix="FlowField",dir=dir, itera=itt,FD=fd)
    ! use keyword argument due to huge number of optional dummy argument
    ! so we need keyword argument to specify the location of actual argument
      call VTR_write_mesh(fd,TGrid%x(:,1),TGrid%y(1,:))
      call VTR_write_var(fd,"Velocity",TVar%u(1:Isize,1:Jsize),TVar%v(1:Isize,1:Jsize))
      call VTR_write_var(fd,"Pressure",TVar%p(1:Isize,1:Jsize))
      call VTR_write_var(fd,"Levelset",TCell%phi(:,:))
      call VTR_write_var(fd,"SolidLvs",TCell%phiS(:,:))
      call VTR_Write_var(fd,"Vof",TCell%vof(:,:))
      call VTR_Write_var(fd,"SolidVof",TCell%vofS(:,:))
      call VTR_close_file(fd)

    !    call VTK_open_file(PREFIX="projectvtk",FD=fd1)
    !    call VTK_write_mesh(FD=fd1,X=x1,Y=y1)
    !    call VTK_write_var(FD=fd1,Name="Velocity",Vx=u,Vy=v)
    !    call VTK_write_var(FD=fd1,Name="Levelset",FIELD=phi)
    !    call VTK_close_file(FD=fd1)

    end subroutine Print_Result_VTK_2D

    SUBROUTINE ReadOldDataPCell(filename,TCell,TVar,FluxP)
      IMPLICIT NONE
      CHARACTER(LEN=80),INTENT(IN):: filename
      TYPE(TVariables),INTENT(INOUT):: TVar
      TYPE(Cell),INTENT(INOUT):: TCell
      REAL(KIND=dp),INTENT(INOUT),DIMENSION(:,:,:),ALLOCATABLE:: FluxP
      INTEGER(KIND=it4b):: i,j
      REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE:: VarRead
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      allocate(VarRead(Isize,Jsize))
      open(unit=5,file=filename,action='read')
      read(5,*)
      read(5,*)
      read(5,*)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TVar%p(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%phi(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%vof(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%phiS(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%vofS(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((FluxP(i,j,1),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((FluxP(i,j,2),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%nxS(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%nyS(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      close(5)
      deallocate(VarRead)
    END SUBROUTINE ReadOldDataPCell

    SUBROUTINE ReadOldDataVelocityCell(filename,TCell,UVVel)
      IMPLICIT NONE
      CHARACTER(LEN=80),INTENT(IN)                          :: filename
      REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: UVVel
      TYPE(Cell),INTENT(INOUT)                              :: TCell
      INTEGER(KIND=it4b)                                    :: i,j
      REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE              :: VarRead
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      allocate(VarRead(Isize,Jsize))
      open(unit=5,file=filename,action='read')
      read(5,*)
      read(5,*)
      read(5,*)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((UVVel(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%phi(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%vof(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%phiS(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((TCell%vofS(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      read(5,"(f24.14)")((VarRead(i,j),i=ibeg,Isize),j=jbeg,Jsize)
      close(5)
    END SUBROUTINE ReadOldDataVelocityCell

    SUBROUTINE ReadOldDataParticle(filename,TraPar,TVar, Lref)
      IMPLICIT NONE
      CHARACTER(LEN=80),INTENT(IN):: filename
      TYPE(TParticle),INTENT(INOUT):: TraPar
      TYPE(TVariables)             :: TVar
      REAL(dp),INTENT(in)          :: Lref
      INTEGER(KIND=it4b)          :: i,j
      CHARACTER(LEN=20)           :: read1,read2,read3,read4,read5
      open(unit=5,file=filename,action='read')
      read(5,*)
      read(5,*) read1,read2,read3,TraPar%np,read5,j
      read(5,"(f24.14)")(TraPar%Posp(i)%x,i=1,TraPar%np)
      read(5,"(f24.14)")(TraPar%Posp(i)%y,i=1,TraPar%np)
      read(5,"(f24.14)")(TraPar%dp(i),i=1,TraPar%np)
      read(5,"(f24.14)")(TraPar%uvp(i)%u,i=1,TraPar%np)
      read(5,"(f24.14)")(TraPar%uvp(i)%v,i=1,TraPar%np)
      close (5)
      do i=1,TraPar%np
        TraPar%Posp(i)%x=TraPar%Posp(i)%x*Lref
        TraPar%Posp(i)%y=TraPar%Posp(i)%y*Lref
        TraPar%dp(i)=TraPar%dp(i)*Lref
        TraPar%uvp(i)%u=TraPar%uvp(i)%u*TVar%Uref
        TraPar%uvp(i)%v=TraPar%uvp(i)%v*TVar%Uref
      end do
    END SUBROUTINE ReadOldDataParticle

    SUBROUTINE setDir(dir_)
      !
      CHARACTER*70,INTENT(in)   :: dir_
      dir = dir_
    END SUBROUTINE setDir

    SUBROUTINE getDir(dir_)
      !
      CHARACTER*70,INTENT(out)   :: dir_
      dir_ = dir
    END SUBROUTINE getDir
End module PrintResult
