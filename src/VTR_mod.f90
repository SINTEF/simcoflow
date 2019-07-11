!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VTR_mod.f90 --- XML VTR ASCII data file
!!
!! Auteur          : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
!! Créé le         : Wed Jul 26 14:36:52 2006
!! Dern. mod. par  : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
!! Dern. mod. le   : Wed Sep 16 14:36:29 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vtr
  use PrecisionVar
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER          :: s=selected_real_kind(6)
  INTEGER, PARAMETER          :: d=selected_real_kind(12)
  character(len=1), PARAMETER :: newline=achar(10)
  INTEGER                     :: iproc=0, nb_procs=1

  TYPE, public :: VTR_file_handle
   PRIVATE
   character(len=80) :: prefix
   INTEGER           :: unit
   INTEGER           :: ni, nj, nk
   INTEGER           :: counter=0
   INTEGER           :: restart=0
   logical           :: first=.true.
  end TYPE VTR_file_handle

  PRIVATE :: handle_error, handle_warning, handle_info

  PRIVATE :: VTR_write_mesh_2d, VTR_write_mesh_3d, &
             VTR_write_vector_2d, VTR_write_vector_3d, &
             VTR_write_scalar_2d, VTR_write_scalar_3d

  public :: VTR_open_file,  &
            VTR_write_mesh, &
            VTR_write_var,  &
            VTR_close_file, &
            VTR_collect_file

  interface VTR_write_mesh
    module procedure VTR_write_mesh_2d, VTR_write_mesh_3d
  end interface

  interface VTR_write_var
    module procedure VTR_write_scalar_2d, VTR_write_scalar_3d, &
                     VTR_write_vector_2d, VTR_write_vector_3d
  end interface

  contains

  subroutine handle_error(name, message)
    IMPLICIT NONE

    character(len=*), INTENT(IN) :: name, message
    INTEGER :: code, errcode=92

    print '(/,"   *** Error *** ", A,": ", A,/)',name, message
    stop
  end subroutine handle_error

  subroutine handle_warning(name, message)
    IMPLICIT NONE
    character(len=*), INTENT(IN) :: name, message

    print '(/,"   *** Warning *** ",A,": ", A,/)',name, message
  end subroutine handle_warning

  subroutine handle_info(name, message)
    IMPLICIT NONE
    character(len=*), INTENT(IN) :: name, message

    print '(/,"   *** Info *** ",A,": ", A,/)',name, message
  end subroutine handle_info

  subroutine VTR_open_file(prefix, dir, proc_rank, num_procs, restart,itera,fd)
    IMPLICIT NONE
    character(len=*), INTENT(IN)         :: prefix
    character(len=*), INTENT(in)         :: dir
    INTEGER, optional, INTENT(IN)        :: proc_rank, num_procs, restart
    INTEGER(kind=it8b),optional,INTENT(IN)          ::itera
    TYPE(VTR_file_handle), INTENT(INOUT) :: fd
    character(len=10) :: rank, snapshot
    character(len=80) :: f
    INTEGER           :: err, iproc
    logical           :: file_opened

    !... Looking for a none connected logical file unit.
    fd%prefix=trim(prefix)
    fd%unit=99
    inquire(unit=fd%unit, opened=file_opened)
    do while (file_opened .and. fd%unit /= 0)

      fd%unit = fd%unit - 1
      inquire(unit=fd%unit, opened=file_opened)

    end do
    if (fd%unit == 0 .and. file_opened) then

      call handle_error("VTR_open_file","All file units from 0 to 99 are already connected.")
      stop

    else

      if ( present(proc_rank) .and. present(num_procs) ) then

        iproc=proc_rank
        nb_procs=num_procs

      else if ( present(proc_rank) ) then

        call handle_error("VTK_open_file","Both PROC_RANK and NUM_PROCS arguments must be present.")

      else if ( present(num_procs) ) then

        call handle_error("VTK_open_file","Both PROC_RANK and NUM_PROCS arguments must be present.")

      end if
      if ((fd%first) .and. (present(restart))) then
         fd%restart=restart
         fd%counter=restart
         fd%first=.false.
      end if
      fd%counter=fd%counter+1
      if(present(itera)) fd%counter = itera
      write(snapshot,'(i8.8)') fd%counter
      if ( present(proc_rank) ) then

        write(rank,'(i8)') iproc
        if(present(itera)) then
          f=trim(adjustl(dir))//'Paraview/'//trim(fd%prefix)//                &
            trim(adjustl(rank))//"_"//trim(adjustl(snapshot))//".vtr"
        end if

      else

        if(present(itera)) then
          f=trim(adjustl(dir))//'Paraview/'//trim(fd%prefix)//"_"//           &
            trim(adjustl(snapshot))//".vtr"
        end if

      end if
      open(unit=fd%unit, file=trim(adjustl(f)), form="FORMATTED", status="replace", &
           action="write", iostat=err)
      if(err /= 0) print '("Problem creating file ",a,".")', trim(f)

    end if
    write(unit=fd%unit, fmt='(100A)') '<VTKFile type="RectilinearGrid" version="0.1" format="ascii">'
 end subroutine VTR_open_file

 subroutine VTR_write_mesh_2d(fd, x, y)
    IMPLICIT NONE

    TYPE(VTR_file_handle), INTENT(INOUT)            :: fd
    REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: x, y
    character(len=10) :: buf1, buf2

    fd%ni=size(x) ; fd%nj=size(y)
    write(buf1,'(I6)') fd%ni
    write(buf2,'(I6)') fd%nj

    write(unit=fd%unit, fmt='(100A)') ' <RectilinearGrid WholeExtent="1 '//trim(adjustl(buf1))//&
                                &' 1 '//trim(adjustl(buf2))//' 1 1 '//'">'
    write(unit=fd%unit, fmt='(100A)') '  <Piece Extent="1 '//trim(adjustl(buf1))//' 1 '//&
                                &trim(adjustl(buf2))//' 1 1 '//'">'
    write(unit=fd%unit, fmt='(100A)') '   <Coordinates>'
    write(unit=fd%unit, fmt='(100A)') '    <DataArray type="Float32" Name="X_COORDINATES" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit, fmt=*) real(x(1:fd%ni),kind=s)
    write(unit=fd%unit, fmt='(100A)') '    </DataArray>'
    write(unit=fd%unit, fmt='(100A)') '    <DataArray type="Float32" Name="Y_COORDINATES" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit, fmt=*) real(y(1:fd%nj),kind=s)
    write(unit=fd%unit, fmt='(100A)') '    </DataArray>'
    write(unit=fd%unit, fmt='(100A)') '    <DataArray type="Float32" Name="Z_COORDINATES" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit, fmt=*) real(0.0_d,kind=s)
    write(unit=fd%unit, fmt='(100A)') '    </DataArray>'
    write(unit=fd%unit, fmt='(100A)') '   </Coordinates>'
    write(unit=fd%unit, fmt='(100A)') '   <PointData>'
  end subroutine VTR_write_mesh_2d

  subroutine VTR_write_mesh_3d(fd, x, y, z)
    IMPLICIT NONE

    TYPE(VTR_file_handle), INTENT(INOUT)            :: fd
    REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: x, y, z
    character(len=10) :: buf1, buf2, buf3

    fd%ni=size(x) ; fd%nj=size(y) ; fd%nk=size(z)
    write(buf1,'(I6)') fd%ni
    write(buf2,'(I6)') fd%nj
    write(buf3,'(I6)') fd%nk

    write(unit=fd%unit, fmt='(100A)') ' <RectilinearGrid WholeExtent="1 '//trim(adjustl(buf1))//&
                                &' 1 '//trim(adjustl(buf2))//' 1 '//trim(adjustl(buf3))//'">'
    write(unit=fd%unit, fmt='(100A)') '  <Piece Extent="1 '//trim(adjustl(buf1))//' 1 '//&
                                &trim(adjustl(buf2))//' 1 '//trim(adjustl(buf3))//'">'
    write(unit=fd%unit, fmt='(100A)') '   <Coordinates>'
    write(unit=fd%unit, fmt='(100A)') '    <DataArray type="Float32" Name="X_COORDINATES" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit, fmt=*) real(x(1:fd%ni),kind=s)
    write(unit=fd%unit, fmt='(100A)') '    </DataArray>'
    write(unit=fd%unit, fmt='(100A)') '    <DataArray type="Float32" Name="Y_COORDINATES" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit, fmt=*) real(y(1:fd%nj),kind=s)
    write(unit=fd%unit, fmt='(100A)') '    </DataArray>'
    write(unit=fd%unit, fmt='(100A)') '    <DataArray type="Float32" Name="Z_COORDINATES" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit, fmt=*) real(z(1:fd%nk),kind=s)
    write(unit=fd%unit, fmt='(100A)') '    </DataArray>'
    write(unit=fd%unit, fmt='(100A)') '   </Coordinates>'
    write(unit=fd%unit, fmt='(100A)') '   <PointData>'
  end subroutine VTR_write_mesh_3d

  subroutine VTR_write_vector_2d(fd, name, vx, vy)
    IMPLICIT NONE

    TYPE(VTR_file_handle), INTENT(IN)                 :: fd
    character(len=*), INTENT(IN)             :: name
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:) :: vx, vy
    REAL(KIND=dp), allocatable, DIMENSION(:,:,:) :: velocity
    INTEGER            :: i, j, code
    character(len=256) :: uname, vname

    if ((size(vx,dim=1) /= fd%ni) .or. &
        (size(vx,dim=2) /= fd%nj) ) call handle_warning("VTR_write_var","Incompatible Vx component and mesh sizes.")

    if ((size(vy,dim=1) /= fd%ni) .or. &
        (size(vy,dim=2) /= fd%nj) ) call handle_warning("VTR_write_var","Incompatible Vy component and mesh sizes.")

    if (.not.Allocated(velocity)) then
      Allocate(velocity(3,1:fd%ni,1:fd%nj),STAT=code)
      if ( code /= 0 ) &
       call handle_error("VTR_write_var","Not enough memory to Allocate VELOCITY array")
    end if

    do j=1, fd%nj
      do i=1, fd%ni

          velocity(1,i,j) = vx(i,j)
          velocity(2,i,j) = vy(i,j)
          velocity(3,i,j) = 0.0_d

      end do
    end do

    write(unit=fd%unit,fmt='(100A)') '        <DataArray type="Float32" Name="'//&
                               &trim(adjustl(name))//'" NumberOfComponents="3" format="ascii">'
    write(unit=fd%unit,fmt=*) real(Velocity(:,1:fd%ni,1:fd%nj),kind=s)
    write(unit=fd%unit,fmt='(100A)') '        </DataArray>'

    uname="X_"//name
    write(unit=fd%unit,fmt='(100A)') '        <DataArray type="Float32" Name="'//&
                               &trim(adjustl(uname))//'" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit,fmt=*) real(velocity(1,1:fd%ni,1:fd%nj),kind=s)
    write(unit=fd%unit,fmt='(100A)') '        </DataArray>'

    vname="Y_"//name
    write(unit=fd%unit,fmt='(100A)') '        <DataArray type="Float32" Name="'//&
                               &trim(adjustl(vname))//'" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit,fmt=*) real(velocity(2,1:fd%ni,1:fd%nj),kind=s)
    write(unit=fd%unit,fmt='(100A)') '        </DataArray>'

    if (Allocated(velocity)) deAllocate(velocity)
  end subroutine VTR_write_vector_2d

  subroutine VTR_write_vector_3d(fd, name, vx, vy, vz)
    IMPLICIT NONE

    TYPE(VTR_file_handle), INTENT(IN)                   :: fd
    character(len=*), INTENT(IN)               :: name
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:,:) :: vx, vy, vz
    REAL(KIND=dp), allocatable, DIMENSION(:,:,:,:) :: velocity
    INTEGER            :: i, j, k, code
    character(len=256) :: uname, vname, wname

    if ((size(vx,dim=1) /= fd%ni) .or. &
        (size(vx,dim=2) /= fd%nj) .or. &
        (size(vx,dim=3) /= fd%nk)) call handle_warning("VTR_write_var","Incompatible Vx component and mesh sizes.")

    if ((size(vy,dim=1) /= fd%ni) .or. &
        (size(vy,dim=2) /= fd%nj) .or. &
        (size(vy,dim=3) /= fd%nk)) call handle_warning("VTR_write_var","Incompatible Vy component and mesh sizes.")

    if ((size(vz,dim=1) /= fd%ni) .or. &
        (size(vz,dim=2) /= fd%nj) .or. &
        (size(vz,dim=3) /= fd%nk)) call handle_warning("VTR_write_var","Incompatible Vz component and mesh sizes.")

    if (.not.Allocated(velocity)) then
      Allocate(velocity(3,1:fd%ni,1:fd%nj,1:fd%nk),STAT=code)
      if ( code /= 0 ) &
       call handle_error("VTR_write_var","Not enough memory to Allocate VELOCITY array")
    end if

    do k=1, fd%nk
      do j=1, fd%nj
        do i=1, fd%ni

          velocity(1,i,j,k) = vx(i,j,k)
          velocity(2,i,j,k) = vy(i,j,k)
          velocity(3,i,j,k) = vz(i,j,k)

        end do
      end do
    end do

    write(unit=fd%unit,fmt='(100A)') '        <DataArray type="Float32" Name="'//&
                               &trim(adjustl(name))//'" NumberOfComponents="3" format="ascii">'
    write(unit=fd%unit,fmt=*) real(Velocity(:,1:fd%ni,1:fd%nj,1:fd%nk),kind=s)
    write(unit=fd%unit,fmt='(100A)') '        </DataArray>'

    uname="X_"//name
    write(unit=fd%unit,fmt='(100A)') '        <DataArray type="Float32" Name="'//&
                               &trim(adjustl(uname))//'" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit,fmt=*) real(velocity(1,1:fd%ni,1:fd%nj,1:fd%nk),kind=s)
    write(unit=fd%unit,fmt='(100A)') '        </DataArray>'

    vname="Y_"//name
    write(unit=fd%unit,fmt='(100A)') '        <DataArray type="Float32" Name="'//&
                               &trim(adjustl(vname))//'" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit,fmt=*) real(velocity(2,1:fd%ni,1:fd%nj,1:fd%nk),kind=s)
    write(unit=fd%unit,fmt='(100A)') '        </DataArray>'

    wname="Z_"//name
    write(unit=fd%unit,fmt='(100A)') '        <DataArray type="Float32" Name="'//&
                               &trim(adjustl(wname))//'" NumberOfComponents="1" format="ascii">'
    write(unit=fd%unit,fmt=*) real(velocity(3,1:fd%ni,1:fd%nj,1:fd%nk),kind=s)
    write(unit=fd%unit,fmt='(100A)') '        </DataArray>'

    if (Allocated(velocity)) deAllocate(velocity)
  end subroutine VTR_write_vector_3d

  subroutine VTR_write_scalar_2d(fd, name, field)
    IMPLICIT NONE

    TYPE(VTR_file_handle), INTENT(IN)                 :: fd
    character(len=*), INTENT(IN)             :: name
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:) :: field

    if ((size(field,dim=1) /= fd%ni) .or. &
        (size(field,dim=2) /= fd%nj)) call handle_warning("VTR_write_var","Incompatible FIELD and MESH sizes.")

      write(unit=fd%unit,fmt='(100A)') '         <DataArray type="Float32" Name="'//&
                                 &trim(adjustl(name))//'" NumberOfComponents="1" format="ascii">'
      write(unit=fd%unit,fmt=*) real(field(1:fd%ni,1:fd%nj),kind=s)
      write(unit=fd%unit,fmt='(100A)') '         </DataArray>'
  end subroutine VTR_write_scalar_2d

  subroutine VTR_write_scalar_3d(fd, name, field)
    IMPLICIT NONE

    TYPE(VTR_file_handle), INTENT(IN)                   :: fd
    character(len=*), INTENT(IN)               :: name
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:,:) :: field

    if ((size(field,dim=1) /= fd%ni) .or. &
        (size(field,dim=2) /= fd%nj) .or. &
        (size(field,dim=3) /= fd%nk)) call handle_warning("VTR_write_var","Incompatible FIELD and MESH sizes.")

      write(unit=fd%unit,fmt='(100A)') '         <DataArray type="Float32" Name="'//&
                                 &trim(adjustl(name))//'" NumberOfComponents="1" format="ascii">'
      write(unit=fd%unit,fmt=*) real(field(1:fd%ni,1:fd%nj,1:fd%nk),kind=s)
      write(unit=fd%unit,fmt='(100A)') '         </DataArray>'
  end subroutine VTR_write_scalar_3d

  subroutine VTR_close_file(fd)
   IMPLICIT NONE

   TYPE(VTR_file_handle), INTENT(INOUT) :: fd
   logical                     :: file_opened

   inquire(unit=fd%unit, opened=file_opened)
   if (file_opened) then
     write(unit=fd%unit,fmt='(100A)') '     </PointData>'
     write(unit=fd%unit,fmt='(100A)') '   </Piece>'
     write(unit=fd%unit,fmt='(100A)') '  </RectilinearGrid>'
     write(unit=fd%unit,fmt='(100A)') '</VTKFile>'
     close(unit=fd%unit)
   else
     call handle_warning("VTR_close_file","No such file to close. Please, check file descriptor.")
   end if
  end subroutine VTR_close_file

  subroutine VTR_collect_file(fd)
    IMPLICIT NONE

    TYPE(VTR_file_handle), INTENT(INOUT) :: fd
    character(len=10) :: rank, snapshot
    character(len=80) :: f, vtrfile
    INTEGER           :: shot, code, err, nt, np, k
    logical           :: file_opened

    !... Looking for a none connected logical file unit.
    if (iproc == 0) then
      fd%unit=99
      inquire(unit=fd%unit, opened=file_opened)
      do while (file_opened .and. fd%unit /= 0)
        fd%unit = fd%unit - 1
        inquire(unit=fd%unit, opened=file_opened)
      end do
      if (fd%unit == 0 .and. file_opened) then
        call handle_warning("VTR_open_file","warning, all file units from 0 to 99 are already connected.")
      else
        f=trim(adjustl(fd%prefix))//".pvd"
        open(unit=fd%unit, file=trim(adjustl(f)), form="FORMATTED", status="replace", &
           action="write", iostat=err)
        if(err /= 0) print '("VTR_collect_file: Error, problem creating file ",a,".")', trim(f)
        write(unit=fd%unit,fmt='(100A)')   '<?xml version="1.0"?>'
        write(unit=fd%unit,fmt='(100A)')   '<VTKFile type="Collection" version="0.1" format="ascii">'
        write(unit=fd%unit,fmt='(100A)')   '  <Collection>'
        nt=len_trim(fd%prefix)
        np=scan(STRING=fd%prefix, SET="/", BACK=.true.)
        vtrfile=fd%prefix(np+1:nt)
        if ( nb_procs == 1 ) then

          do shot = 1, fd%counter
            write(snapshot, '(i6)') shot
            write(unit=fd%unit,fmt='(100A)') '    <DataSet timestep="'//trim(adjustl(snapshot))//&
                                      &'" part="0'//'" file="'//trim(adjustl(vtrfile))//&
                                      &"_"//trim(adjustl(snapshot))//'.vtr"/>'
          end do

        else

          do k = 0, nb_procs-1
            write(rank, '(i6)') k
            do shot = 1, fd%counter
              write(snapshot, '(i6)') shot
              write(unit=fd%unit,fmt='(100A)') '    <DataSet timestep="'//trim(adjustl(snapshot))//&
                                        &'" part="'//trim(adjustl(rank))//'" file="'//&
                                        &trim(adjustl(vtrfile))//"_"//trim(adjustl(rank))//&
                                        &"_"//trim(adjustl(snapshot))//'.vtr"/>'
            end do
          end do

        end if
        write(unit=fd%unit,fmt='(100A)')    '  </Collection>'
        write(unit=fd%unit,fmt='(100A)')    '</VTKFile>'
        close(unit=fd%unit)
      end if
    end if
    fd%counter=0 ; fd%restart=0; fd%first=.true. ; iproc=0 ; nb_procs=1
  end subroutine VTR_collect_file
end module vtr
