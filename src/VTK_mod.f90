!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VTK_mod.f90 --- VTK data file format
!!
!! Auteur          : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
!! Cr�� le         : Wed Jul 26 14:36:52 2006
!! Dern. mod. par  : Jalel Chergui (LIMSI-CNRS) <Jalel.Chergui@limsi.fr>
!! Dern. mod. le   : Wed Sep 16 14:35:06 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vtk
  use PrecisionVar
  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER          :: s=selected_real_kind(6)
  INTEGER, PARAMETER          :: d=selected_real_kind(12)
  character(len=1), PARAMETER :: newline=achar(10)
  INTEGER                     :: iproc=0, nb_procs=1
  TYPE, PUBLIC :: VTK_file_handle
    PRIVATE
    character(len=80) :: prefix
    INTEGER           :: unit
    INTEGER           :: ni, nj, nk
    INTEGER           :: counter=0
    INTEGER           :: restart=0
    logical           :: first=.true.
  end TYPE VTK_file_handle

  PRIVATE :: handle_error, handle_warning, handle_info

  PRIVATE :: VTK_write_mesh_2d,   VTK_write_mesh_3d,   &
             VTK_write_vector_2d, VTK_write_vector_3d, &
             VTK_write_scalar_2d, VTK_write_scalar_3d

  PUBLIC :: VTK_open_file,  &
            VTK_write_mesh, &
            VTK_write_var,  &
            VTK_close_file, &
            VTK_collect_file

  interface VTK_write_mesh
    module procedure VTK_write_mesh_2d, VTK_write_mesh_3d
  end interface

  interface VTK_write_var
    module procedure VTK_write_scalar_2d, VTK_write_scalar_3d, &
                     VTK_write_vector_2d, VTK_write_vector_3d
  end interface

  contains

  subroutine handle_error(name, message)
    IMPLICIT NONE

    character(len=*), INTENT(IN) :: name, message

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

  subroutine VTK_open_file(prefix, proc_rank, num_procs, restart, fd)
    IMPLICIT NONE

    character(len=*), INTENT(IN)         :: prefix
    INTEGER, optional, INTENT(IN)        :: proc_rank, num_procs, restart
    TYPE(VTK_file_handle), INTENT(INOUT) :: fd
    character(len=10) :: rank, snapshot
    character(len=80) :: f
    character(len=256):: MAIN_header
    INTEGER           :: err
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

      call handle_error("VTK_open_file","All file units from 0 to 99 are already connected.")
      stop

    else

      if ( present(proc_rank) .and. present(num_procs) ) then

        iproc    = proc_rank
        nb_procs = num_procs

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
      write(snapshot,'(i8.8)') fd%counter
      if ( present(proc_rank) ) then
        write(rank,'(i8)') iproc
        f=trim(fd%prefix)//"_"//trim(adjustl(rank))//"_"//trim(adjustl(snapshot))//".vtk"
      else
        f=trim(fd%prefix)//"_"//trim(adjustl(snapshot))//".vtk"
      end if
      open(unit=fd%unit, file=trim(adjustl(f)), form="UNFORMATTED", access="STREAM", status="replace", &
           action="write", iostat=err)
      if(err /= 0) print '("Problem creating file ",a,".")', trim(f)
    end if
    MAIN_header="# vtk DataFile Version 2.0"//newline//"(c) J.C. April 2009"//newline//"BINARY"//newline
    write(unit=fd%unit) trim(MAIN_header)
  end subroutine VTK_open_file

  subroutine VTK_write_mesh_2d(fd, x, y)
    IMPLICIT NONE

    TYPE(VTK_file_handle), INTENT(INOUT)   :: fd
    REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: x, y
    character(len=30)  :: buf1, buf2
    character(len=256) :: GRID_header
    INTEGER, PARAMETER :: nk=1

    fd%ni=size(x) ; fd%nj=size(y)
    write(buf1,'(i8," ",i8," ",i8)') fd%ni,fd%nj,fd%nk
    GRID_header="DATASET RECTILINEAR_GRID"//newline//"DIMENSIONS "//trim(adjustl(buf1))//newline
    write(unit=fd%unit) trim(GRID_header)
    write(buf2,'(i8)') fd%ni
    GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(x(1:fd%ni),kind=s),newline
    write(buf2,'(i8)') fd%nj
    GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(y(1:fd%nj),kind=s),newline
    write(buf2,'(i8)') fd%nk
    GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(0.0_d,kind=s),newline
    write(buf2,'(i8)') fd%ni*fd%nj
    GRID_header="POINT_DATA "//trim(adjustl(buf2))//newline
    write(unit=fd%unit) trim(GRID_header)
  end subroutine VTK_write_mesh_2d

  subroutine VTK_write_mesh_3d(fd, x, y, z)
    IMPLICIT NONE

    TYPE(VTK_file_handle), INTENT(INOUT)   :: fd
    REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: x, y, z
    character(len=30)  :: buf1, buf2
    character(len=256) :: GRID_header

    fd%ni=size(x) ; fd%nj=size(y) ; fd%nk=size(z)
    write(buf1,'(i8," ",i8," ",i8)') fd%ni,fd%nj,fd%nk
    GRID_header="DATASET RECTILINEAR_GRID"//newline//"DIMENSIONS "//trim(adjustl(buf1))//newline
    write(unit=fd%unit) trim(GRID_header)
    write(buf2,'(i8)') fd%ni
    GRID_header="X_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(x(1:fd%ni),kind=s),newline
    write(buf2,'(i8)') fd%nj
    GRID_header="Y_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(y(1:fd%nj),kind=s),newline
    write(buf2,'(i8)') fd%nk
    GRID_header="Z_COORDINATES "//trim(adjustl(buf2))//" float"//newline
    write(unit=fd%unit) trim(GRID_header),real(z(1:fd%nk),kind=s),newline
    write(buf2,'(i8)') fd%ni*fd%nj*fd%nk
    GRID_header="POINT_DATA "//trim(adjustl(buf2))//newline
    write(unit=fd%unit) trim(GRID_header)
  end subroutine VTK_write_mesh_3d

  subroutine VTK_write_vector_2d(fd, name, vx, vy)
    IMPLICIT NONE

    TYPE(VTK_file_handle), INTENT(IN)        :: fd
    character(len=*), INTENT(IN)             :: name
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:) :: vx, vy
    REAL(KIND=dp), allocatable, DIMENSION(:,:,:) :: velocity
    INTEGER                                     :: i, j, code=0
    character(len=256)                          :: uname, vname, VAR_header

    if ((size(vx,dim=1) /= fd%ni) .or. &
        (size(vx,dim=2) /= fd%nj)) call handle_warning("VTK_write_var","Incompatible X component and mesh sizes.")

    if ((size(vy,dim=1) /= fd%ni) .or. &
        (size(vy,dim=2) /= fd%nj)) call handle_warning("VTK_write_var","Incompatible Y component and mesh sizes.")

    if (.not.Allocated(velocity)) then

      Allocate(velocity(3,fd%ni,fd%nj),STAT=code)
      if ( code /= 0 ) &
      call handle_error("VTK_write_var","Not enough memory to Allocate VELOCITY array")

    end if

    do j=1, fd%nj
      do i=1, fd%ni

          velocity(1,i,j) = vx(i,j)
          velocity(2,i,j) = vy(i,j)
          velocity(3,i,j) = 0.0_d

      end do
    end do

    VAR_header="VECTORS "//trim(adjustl(name))//" float "//newline
    write(unit=fd%unit) trim(VAR_header),real(velocity(:,1:fd%ni,1:fd%nj),kind=s),newline

    uname="X_"//name
    VAR_header="SCALARS "//trim(adjustl(uname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(1,1:fd%ni,1:fd%nj),kind=s),newline

    vname="Y_"//name
    VAR_header="SCALARS "//trim(adjustl(vname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(2,1:fd%ni,1:fd%nj),kind=s),newline

    if (Allocated(velocity)) deAllocate(velocity)
  end subroutine VTK_write_vector_2d

  subroutine VTK_write_vector_3d(fd, name, vx, vy, vz)
    IMPLICIT NONE

    TYPE(VTK_file_handle), INTENT(IN)          :: fd
    character(len=*), INTENT(IN)               :: name
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:,:) :: vx, vy, vz
    REAL(KIND=dp), allocatable, DIMENSION(:,:,:,:) :: velocity
    INTEGER                                       :: i, j, k, code=0
    character(len=256)                            :: uname, vname, wname, VAR_header

    if ((size(vx,dim=1) /= fd%ni) .or. &
        (size(vx,dim=2) /= fd%nj) .or. &
        (size(vx,dim=3) /= fd%nk)) call handle_warning("VTK_write_var","Incompatible X component and mesh sizes.")

    if ((size(vy,dim=1) /= fd%ni) .or. &
        (size(vy,dim=2) /= fd%nj) .or. &
        (size(vy,dim=3) /= fd%nk)) call handle_warning("VTK_write_var","Incompatible Y component and mesh sizes.")

    if ((size(vz,dim=1) /= fd%ni) .or. &
        (size(vz,dim=2) /= fd%nj) .or. &
        (size(vz,dim=3) /= fd%nk)) call handle_warning("VTK_write_var","Incompatible Z component and mesh sizes.")

    if (.not.Allocated(velocity)) then
      Allocate(velocity(3,fd%ni,fd%nj,fd%nk),STAT=code)
      if ( code /= 0 ) &
      call handle_error("VTK_write_var","Not enough memory to Allocate VELOCITY array")
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

    VAR_header="VECTORS "//trim(adjustl(name))//" float "//newline
    write(unit=fd%unit) trim(VAR_header),real(velocity(:,1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline

    uname="X_"//name
    VAR_header="SCALARS "//trim(adjustl(uname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(1,1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline

    vname="Y_"//name
    VAR_header="SCALARS "//trim(adjustl(vname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(2,1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline

    wname="Z_"//name
    VAR_header="SCALARS "//trim(adjustl(wname))//" float "//newline//"LOOKUP_TABLE default"//newline
      write(unit=fd%unit) trim(VAR_header),real(velocity(3,1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline

    if (Allocated(velocity)) deAllocate(velocity)
  end subroutine VTK_write_vector_3d

  subroutine VTK_write_scalar_2d(fd, name, field)
    IMPLICIT NONE

    TYPE(VTK_file_handle), INTENT(IN)        :: fd
    character(len=*), INTENT(IN)             :: name
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:) :: field
    character(len=256) :: VAR_header

    if ((size(field,dim=1) /= fd%ni) .or. (size(field,dim=2) /= fd%nj)) &
       call handle_warning("VTK_write_var","Incompatible FIELD and MESH sizes.")

    VAR_header="SCALARS "//trim(adjustl(name))//" float "//newline//"LOOKUP_TABLE default"//newline
    write(unit=fd%unit) trim(VAR_header),real(field(1:fd%ni,1:fd%nj),kind=s),newline
  end subroutine VTK_write_scalar_2d

  subroutine VTK_write_scalar_3d(fd, name, field)
    IMPLICIT NONE

    TYPE(VTK_file_handle), INTENT(IN)          :: fd
    character(len=*), INTENT(IN)               :: name
    REAL(KIND=dp), INTENT(IN), DIMENSION(:,:,:) :: field
    character(len=256) :: VAR_header

    if ((size(field,dim=1) /= fd%ni) .or. &
        (size(field,dim=2) /= fd%nj) .or. &
        (size(field,dim=3) /= fd%nk)) call handle_warning("VTK_write_var","Incompatible FIELD and MESH sizes.")

    VAR_header="SCALARS "//trim(adjustl(name))//" float "//newline//"LOOKUP_TABLE default"//newline
    write(unit=fd%unit) trim(VAR_header),real(field(1:fd%ni,1:fd%nj,1:fd%nk),kind=s),newline
  end subroutine VTK_write_scalar_3d

  subroutine VTK_close_file(fd)
   IMPLICIT NONE

   TYPE(VTK_file_handle), INTENT(IN) :: fd
   logical :: file_opened

   inquire(unit=fd%unit, opened=file_opened)
   if (file_opened) then
     close(unit=fd%unit)
   else
     call handle_warning("VTK_close_file","No such file to close. Please, check file descriptor.")
   end if
  end subroutine VTK_close_file

  subroutine VTK_collect_file(fd)
    IMPLICIT NONE

    TYPE(VTK_file_handle), INTENT(INOUT) :: fd
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
        call handle_warning("VTK_open_file","warning, all file units from 0 to 99 are already connected.")
      else
        f=trim(adjustl(fd%prefix))//".pvd"
        open(unit=fd%unit, file=trim(adjustl(f)), form="FORMATTED", status="replace", &
           action="write", iostat=err)
        if(err /= 0) print '("VTK_collect_file: Error, problem creating file ",a,".")', trim(f)
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
    fd%counter=0 ; fd%restart=0 ; fd%first=.true. ; iproc=0 ; nb_procs=1
  end subroutine VTK_collect_file
end module vtk
