Module Solver
    USE PrecisionVar
    USE Mesh, ONLY : TsimcoMesh, getMeshSizes, Grid, Cell
    USE StateVariables, ONLY : TVariables, TWave, getSolverVariables, Boundary_Condition_Var
    USE Constants, ONLY : g, epsi
    USE CutCell, ONLY : Grid_Preprocess, NewCellFace
    USE Clsvof, ONLY : SolidObject, Coupled_LS_VOF, Initial_ClsVofUV, SolidVolumeFraction, ComputeForceObject
    USE PrintResult, ONLY : getDir, ReadOldDataPCell, ReadOldDataVelocityCell, ReadOldDataParticle
    USE PrintResult, ONLY : Print_Result_Tecplot_UCent, Print_Result_Tecplot_VCent, Print_Result_VTK_2D
    USE PrintResult, ONLY : Print_Result_Tecplot_PCent
    USE ComputePUV, ONLY : UpdatePUV
    USE MPI
    USE Particles, ONLY : TParticle
    IMPLICIT NONE
    PRIVATE
    INTEGER(it8b) :: IttBegin
    TYPE,PUBLIC:: SolverTime
      INTEGER(kind=it8b):: iter
      REAL(KIND=dp):: cfl
      REAL(KIND=dp):: dt,PhysT,NondiT
    End TYPE SolverTime
    TYPE,PUBLIC:: SolverConvergence
      REAL(KIND=dp):: N1,N2,NInf,N1c,N2c,NInfc
    End TYPE SolverConvergence
    REAL(KIND=dp),PARAMETER:: Twalls = 300.d0
    REAL(KIND=dp),PARAMETER:: Twall = 400.d0
    REAL(KIND=dp),PRIVATE:: tprint(7),xprintex(0:7)
    INTEGER(kind=it4b),PRIVATE:: WavePrint
    PUBLIC:: IterationSolution
    Interface IterationSolution
      Module Procedure IterationSolution
    End Interface IterationSolution
    Contains
    subroutine IterationSolution(simcomesh, TVar, wave,    &
                                 TraPar,BoomCase,iprint)
      IMPLICIT NONE
      TYPE(TsimcoMesh) , INTENT(inout) :: simcomesh
      TYPE(TVariables),INTENT(INOUT):: TVar
      TYPE(TWave), INTENT(in) :: wave
      INTEGER(kind=it4b),INTENT(IN):: iprint
      TYPE(TParticle),INTENT(INOUT):: TraPar
      TYPE(SolidObject),INTENT(INOUT):: BoomCase
      TYPE(Cell):: PCellO,UCellO,VCellO
      TYPE(TsimcoMesh) :: simcomesh0
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: GraP
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: Flux_n1
      TYPE(SolverTime):: Time
      TYPE(SolverConvergence):: UConv,VConv,PConv
      TYPE(TVariables):: TVar_n
      INTEGER(kind=it4b):: iprint1,i,j
      INTEGER(kind=it8b):: itt,tempvel
      REAL(KIND=dp):: velaver,timeb
      CHARACTER(LEN=80):: filename,curd
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      INTEGER(it8b) :: IttRun
      LOGICAL :: RunAgain
      CHARACTER*70 :: dir
      call getSolverVariables(IttRun_=IttRun, RunAgain_=RunAgain)
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      call getDir(dir)
      allocate(Flux_n1(0:Isize+1,0:Jsize+1,2))
      Flux_n1(:,:,:)=0.d0
      allocate(GraP(1:Isize,1:Jsize))
      Time%iter=10**6
      Time%NondiT=0.d0
      Time%Cfl=0.3d0
      iprint1=iprint
      tempvel=0
      velaver=0.d0
      timeb=20.d0
    ! Set up time for pring wave profile
   !   tprint(1)=0.50074d0
   !   tprint(2)=0.75264d0
   !   tprint(3)=0.8018d0
   !   tprint(4)=0.85095d0
   !   tprint(5)=0.89396d0
   !   tprint(6)=1.02912d0
   !   tprint(7)=1.20116d0
    ! Set up time for Overtaking collisions
      tprint(1)=2.59892d0
      tprint(2)=3.50209d0
      tprint(3)=4.14721d0
      tprint(4)=4.6971d0
      tprint(5)=5.5972d0
      tprint(6)=6.60174d0
      tprint(7)=8.39887d0
    ! Set up printing position for overtaking collisions
      xprintex(0)=0.d0
      xprintex(1)=1.301585d0
      xprintex(2)=1.294936d0
      xprintex(3)=1.3711241d0 !1.462082d0 ! for lowest point between two waves
      xprintex(4)=1.456825d0
      xprintex(5)=1.46058d0
      xprintex(6)=1.469767d0
      xprintex(7)=1.495967d0
    ! Flag to print wave profile (1 for printing)
    ! This variables is assigned value in ComputeTimeStep Subroutine
      WavePrint=0
      call PrintWaterWave(Time%PhysT,simcomesh%PGrid,simcomesh%PCell)
      simcomesh0 = TsimcoMesh(Isize, Jsize)
    ! Calculate threshold for MUSCL limiter
      if(RunAgain.eqv..TRUE.) then
        Write(curd,'(i8.8)') IttRun
        filename=trim(adjustl(dir))//'Tecplot/Pressure_'//trim(curd)//'.dat'
        call ReadOldDataPCell(filename,simcomesh%PCell,TVar,Flux_n1)
        filename=trim(adjustl(dir))//'Tecplot/Uvelocity_'//trim(curd)//'.dat'
        call ReadOldDataVelocityCell(filename,simcomesh%UCell,TVar%u)
        filename=trim(adjustl(dir))//'Tecplot/Vvelocity_'//trim(curd)//'.dat'
        call ReadOldDataVelocityCell(filename,simcomesh%VCell,TVar%v)
        filename=trim(adjustl(dir))//'Tecplot/Particles_'//trim(curd)//'.dat'
        call ReadOldDataParticle(filename,TraPar,TVar, simcomesh%PGrid%Lref)
        filename=trim(adjustl(dir))//'Convergence.dat'
        call ReadFileConvergence(filename,Time,PConv,BoomCase)
        BoomCase%XBar1=BoomCase%Posp%x-BoomCase%Wobj/2.d0
        BoomCase%XBar2=BoomCase%Posp%x+BoomCase%Wobj/2.d0
        BoomCase%LBar=1.5d0/simcomesh%PGrid%Lref
        BoomCase%YBar=BoomCase%Posp%y-dsqrt((BoomCase%Dobj/2.d0)**2.d0-        &
                 (BoomCase%Wobj/2.d0)**2.d0)-BoomCase%LBar

        call SolidVolumeFraction(simcomesh%PGrid,simcomesh%PCell,BoomCase)
        call SolidVolumeFraction(simcomesh%UGrid,simcomesh%UCell,BoomCase)
        call SolidVolumeFraction(simcomesh%VGrid,simcomesh%VCell,BoomCase)

        call Grid_Preprocess(simcomesh,TVar,IttRun)
        call NewCellFace(simcomesh)
 !       call Boundary_Condition_Var(PGrid,TVar,Time%NondiT)
        Time%PhysT = Time%Nondit*simcomesh%PGrid%Lref/TVar%URef
        IttBegin=IttRun+1
!        BoomCase%us=0.d0
!        BoomCase%vs=0.d0
      else
        IttBegin=1
      end if
      do itt = IttBegin,Time%iter
!        if(itt>=610) then
!          print*,itt
!          iprint1=1
!        end if
       ! if(itt==365) pause
        call AdamBasforthCrankNicolson(simcomesh, simcomesh0,  &
                 wave, TVar,UConv,VConv,PConv,Time,Flux_n1,TraPar, &
                 BoomCase,itt)
        if(mod(itt,TraPar%IParInlet)==0) then
          call TraPar%ParticleInletCondition(simcomesh%PGrid,simcomesh%PCell,wave)
        end if
        call TraPar%TrackingParticles(simcomesh%PGrid,simcomesh%PCell,BoomCase,TVar,Time%dt)
        Time%NondiT = Time%NondiT+Time%dt
        Time%PhysT = Time%Nondit*simcomesh%PGrid%Lref/TVar%URef
        if(itt==1) then
          open(unit=5,file=trim(adjustl(dir))//'Convergence.dat')
          close(5,status='delete')
          open(unit=5,file=trim(adjustl(dir))//'WaveJetVelocity.dat')
          close(5,status='delete')
          open(unit=5,file=trim(adjustl(dir))//'WaterHeight.dat')
          close(5,status='delete')
          open(unit=5,file=trim(adjustl(dir))//'WaterPressure.dat')
          close(5,status='delete')
          open(unit=5,file=trim(adjustl(dir))//'ObjectForce.dat')
          close(5,status='delete')
        end if

    !    call PrintWaterJetVelocity(Time%NonDiT,PGrid,PCell,TVar,itt,tempvel,   &
    !                                                              velaver,timeb)
    !    call PrintWaveFront(Time%NonDiT,PGrid,PCell)
    !    call PrintWaterHeight(Time%NonDiT,I1,I2,I3,I4,PGrid,PCell)
    !    call PrintWaterPressure(Time%NonDiT,J1,J2,J3,J4,TVar)
        call PrintHistory(itt,Time,PConv,BoomCase)
        if(mod(itt,iprint1)==0)then
          write(*,*) itt,Time%PhysT,Time%NondiT
          call Print_Result_Tecplot_PCent(simcomesh%PGrid,TVar,simcomesh%PCell,TraPar,Flux_n1,itt,1)
          call Print_Result_Tecplot_UCent(simcomesh%UGrid,TVar,simcomesh%UCell,itt)
          call Print_Result_Tecplot_VCent(simcomesh%VGrid,TVar,simcomesh%VCell,itt)
          call Print_Result_VTK_2D(simcomesh%PGrid,TVar,simcomesh%PCell,itt)
        end if
        if(WavePrint/=0) then
          call PrintWaterWave(Time%PhysT,simcomesh%PGrid,simcomesh%PCell)
          WavePrint=0
        end if
      end do
      deallocate(GraP,Flux_n1)
    end subroutine IterationSolution

    SUBROUTINE AdamBasforthCrankNicolson(simcomesh, simcomesh0, wave, &
               TVar,UConv,VConv,PConv,Time,Flux_n1,  &
               TraPar,BoomCase,itt)
      IMPLICIT NONE
      TYPE(TsimcoMesh), INTENT(inout) :: simcomesh
      TYPE(TsimcoMesh), INTENT(inout) :: simcomesh0
      TYPE(TWave), INTENT(in)         :: wave
      TYPE(TVariables),INTENT(INOUT):: TVar
      TYPE(SolverTime),INTENT(INOUT):: Time
      TYPE(TParticle),INTENT(INOUT):: TraPar
      TYPE(SolverConvergence),INTENT(OUT):: UConv,VConv,PConv
      TYPE(SolidObject),INTENT(INOUT):: BoomCase
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(INOUT):: Flux_n1
      INTEGER(kind=it8b),INTENT(IN):: itt
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: SPar
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: SParU,SParV,VolPar,VolParU,VolParV
      REAL(KIND=dp):: dt,Se,ForceObj
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      INTEGER(it8b) :: IttRun
      LOGICAL :: RunAgain
      CHARACTER*70 :: dir
      call getSolverVariables(IttRun_=IttRun, RunAgain_=RunAgain)
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      call getDir(dir)
      if(itt==1) then
        BoomCase%us=0.d0
        BoomCase%vs=0.d0
      end if
      allocate(SPar(Isize,Jsize,2))
      allocate(SParU(Isize,Jsize))
      allocate(SParV(Isize,Jsize))
      allocate(VolPar(Isize,Jsize))
      allocate(VolParU(Isize,Jsize))
      allocate(VolParV(Isize,Jsize))
      ! First Runge-Kutta substep
      call CopyNewCell(simcomesh0%PCell,simcomesh%PCell)
      call CopyNewCell(simcomesh0%UCell,simcomesh%UCell)
      call CopyNewCell(simcomesh0%VCell,simcomesh%VCell)
      call ComputeForceObject(BoomCase,simcomesh%PGrid,simcomesh%PCell,simcomesh%VCell,TVar,ForceObj)
   !   BoomCase%asy=-(1.14d0*omew)**2.d0*Amp0*dsin(1.14d0*omew*Time%NondiT+pi/2.d0)/    &
   !                 (TVar%Uref**2.d0/Lref)
      if(RunAgain.eqv..FALSE.) then
        BoomCase%asy=(ForceObj-BoomCase%Mobj*g)/BoomCase%Mobj
      elseif(itt>ittRun+1) then
        BoomCase%asy=(ForceObj-BoomCase%Mobj*g)/BoomCase%Mobj
      end if
      BoomCase%asx=0.d0
      if(itt<10) then
        BoomCase%us=0.d0
        BoomCase%vs=0.d0
        BoomCase%asy=0.d0
      end if
      call ComputeTimeStep(simcomesh%UGrid,simcomesh%VGrid,TVar,BoomCase,Time)
      open(unit=10,file=trim(adjustl(dir))//'ObjectForce.dat',access='append')
      write(10,*)Time%NondiT,ForceObj-BoomCase%Mobj*g,BoomCase%asy,itt
      close(10)
      dt=time%dt
      if(itt>1) then
        call Coupled_LS_VOF(simcomesh%PGrid,simcomesh%PCell,simcomesh%UCell,simcomesh%VCell,TVar,wave, BoomCase,       &
                                                             Time%NondiT,dt,itt)
        call Initial_ClsVofUV(simcomesh%PCell,simcomesh%PGrid,simcomesh%UCell,simcomesh%UGrid,VolPar,SPar,VolParU,     &
                                                             SParU,BoomCase,0)
        call Initial_ClsVofUV(simcomesh%PCell,simcomesh%PGrid,simcomesh%VCell,simcomesh%VGrid,VolPar,SPar,VolParV,     &
                                                             SParV,BoomCase,1)
        call Grid_Preprocess(simcomesh,TVar,itt)
        call NewCellFace(simcomesh)
        call Boundary_Condition_Var(simcomesh%PGrid,simcomesh%PCell,TVar,wave, Time%NondiT)
        call InterNewVar(simcomesh0%PCell,simcomesh0%UCell,simcomesh0%VCell,simcomesh%PCell,simcomesh%UCell,&
             &           simcomesh%VCell,simcomesh%PGrid,TVar,    &
                                                                BoomCase%vs)
      else
        SParU(:,:)=0.d0;SParV(:,:)=0.d0
        VolParU(:,:)=0.d0;VolParV(:,:)=0.d0
      end if
      call UpdatePUV(simcomesh%UGrid,simcomesh%VGrid,simcomesh%PGrid,simcomesh0%PCell,simcomesh0%UCell, &
           &         simcomesh0%VCell,simcomesh%PCell,simcomesh%UCell,       &
           simcomesh%VCell,TVar,Flux_n1,TraPar,VolParU,VolParV,SParU,SParV,BoomCase,dt,itt)
      ! Calculate the three kind of norm for convergence
      deallocate(SPar)
      deallocate(SParU,SParV)
      deallocate(VolPar)
      deallocate(VolParU,VolParV)
    END SUBROUTINE AdamBasforthCrankNicolson

    SUBROUTINE ComputeTimeStep(UGrid,VGrid,TVar,BoomCase,Time)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: UGrid,VGrid
      TYPE(TVariables),INTENT(IN):: TVar
      TYPE(SolidObject),INTENT(IN):: BoomCase
      TYPE(SolverTime),INTENT(OUT):: Time
      INTEGER(kind=it4b):: i,j
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      Time%dt = 1.d0
      do j = jbeg,jbeg+Jsize-1
        do i = ibeg,ibeg+Isize-1
          Time%dt=dmin1(Time%dt,Time%cfl*Ugrid%dx(i,j)/dabs(TVar%u(i,j)))
          Time%dt=dmin1(Time%dt,Time%cfl*VGrid%dy(i,j)/dabs(TVar%v(i,j)))
          Time%dt=dmin1(Time%dt,2.d0*Time%cfl*VGrid%dy(i,j)/(dabs(TVar%v(i,j))+&
                         dsqrt(TVar%v(i,j)**2.d0+4.d0*VGrid%dy(i,j)/TVar%Fr)),      &
                                Time%cfl*Vgrid%dy(i,j)/(TVar%Vint/TVar%Uref))
          Time%dt=dmin1(Time%dt,Time%cfl*Vgrid%dy(i,j)/(TVar%Vint/TVar%Uref))
          Time%dt=dmin1(Time%dt,Time%cfl*Vgrid%dy(i,j)/dabs(BoomCase%vs))
          Time%dt=dmin1(Time%dt,(-BoomCase%vs+                                 &
             dsqrt(dmax1(dabs(BoomCase%vs**2.d0+2.d0*dabs(BoomCase%asy)*0.05d0* &
             VGrid%dy(i,j)),1d-30)))/dabs(BoomCase%asy+1.d-20))
        end do
        Time%dt=dmin1(Time%dt,Time%cfl*Ugrid%dx(1,j)/dabs(TVar%u(0,j)))
      end do
    END SUBROUTINE ComputeTimeStep

    SUBROUTINE PrintHistory(itt,Time,TNorm,BoomCase)
      IMPLICIT NONE
      INTEGER(kind=it8b),INTENT(IN)     :: itt
      TYPE(SolverTime),INTENT(IN)       :: Time
      TYPE(SolverConvergence),INTENT(IN):: TNorm
      TYPE(SolidObject),INTENT(IN)      :: BoomCase
      CHARACTER*70 :: dir
      call getDir(dir)
      open(unit=5,file=trim(adjustl(dir))//'Convergence.dat',access='append')
      write(5,76) itt,Time%NondiT,TNorm%N1,TNorm%N2,TNorm%Ninf,TNorm%N1c,      &
                     TNorm%N2c,TNorm%Ninfc,BoomCase%Posp%y,                    &
                     BoomCase%PospO%y,BoomCase%vs,BoomCase%asy
      close(5)
 76	  format(I10,11(f24.14))
    END SUBROUTINE

    SUBROUTINE InterNewVar(PCellO,UCellO,VCellO,PCell,UCell,VCell,PGrid,TVar,vb)
       IMPLICIT NONE
       TYPE(Cell),INTENT(IN):: PCellO,UCellO,VCellO,PCell,UCell,VCell
       TYPE(TVariables),INTENT(INOUT):: TVar
       TYPE(Grid),INTENT(IN):: PGrid
       REAL(KIND=dp),intent(in):: vb
       INTEGER(KIND=it4b):: i,j,ii,jj,temp
       REAL(KIND=dp):: Pu,Pv,SumP,vfa,uleft,uright,vbot,vtop
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
       do i = ibeg,ibeg+Isize-1
         do j = jbeg,jbeg+Jsize-1
           ! For pressure cell
           ! from momentum exchange cell to normal cut cell
           if(UCellO%MoExCell(i,j)==1.and.UCell%MoExCell(i,j)==0) then
             ii = UCellO%MsCe(i,j,1)
             jj = UCellO%MsCe(i,j,2)
            ! TVar%u(ii,jj) = TVar%u(ii,jj)*(UCellO%vof(i,j)+UCellO%vof(ii,jj))/  &
            !                 (UCellO%vof(ii,jj)+UCellO%vof(i,j)*UCell%delh(i,j)/&
            !                 UCell%delh(ii,jj))
            ! TVar%u(ii,jj) = (TVar%u(ii,jj)*UCellO%vof(ii,jj)+TVar%u(i,j)*      &
            !                 UCellO%vof(i,j)*(1.d0-UCell%delh(i,j)/            &
            !                 UCell%delh(ii,jj)))/UCellO%vof(ii,jj)
            !  TVar%u(i,j) = TVar%u(i,j)*UCell%delh(i,j)/UCell%delh(ii,jj)
           end if
           ! from normal cut cell to momentum exchange cell
           if(UCell%MoExCell(i,j)==1.and.UCellO%MoExCell(i,j)==0) then
             ii = UCell%MsCe(i,j,1)
             jj = UCell%MsCe(i,j,2)
            ! TVar%u(ii,jj) = (TVar%u(i,j)*UCellO%vof(i,j)+TVar%u(ii,jj)*        &
            !                UCellO%vof(ii,jj))/(UCellO%vof(i,j)+UCellO%vof(ii,jj))
            ! TVar%u(i,j) = TVar%u(ii,jj)
           end if
           ! for U-velocity (New velocity cell)
           if(UCellO%VofS(i,j)>1.d0-epsi.and.UCell%VofS(i,j)<1.d0-epsi) then
             ii = UCell%MsCe(i,j,1)
             jj = UCell%MsCe(i,j,2)
           ! print*,'new u',i,j
             TVar%u(i,j) = 0.d0 !TVar%u(ii,jj)
             if(PCell%EEdge_Area(i,j)>epsi) then
               uleft=((-TVar%v(i,j)*PCell%NEdge_Area(i,j)+                     &
                        TVar%v(i,j-1)*PCell%SEdge_Area(i,j))*                  &
                        PGrid%dx(i,j)+TVar%u(i-1,j)*                           &
                        PCell%WEdge_Area(i,j)*PGrid%dy(i,j)+                   &
                        vb*PCell%nyS(i,j)*PCell%Wllh(i,j))/                    &
                       (PCell%EEdge_Area(i,j)+1.d-30)/PGrid%dy(i,j)
             else
               uleft=0.d0
             end if
             if(PCell%WEdge_Area(i+1,j)>epsi) then
               uright=((TVar%v(i+1,j)*PCell%NEdge_Area(i+1,j)-                 &
                        TVar%v(i+1,j-1)*PCell%SEdge_Area(i+1,j))*              &
                        PGrid%dx(i+1,j)+TVar%u(i+1,j)*                         &
                        PCell%EEdge_Area(i+1,j)*PGrid%dy(i+1,j)-               &
                        vb*PCell%nyS(i+1,j)*PCell%Wllh(i+1,j))/                &
                        (PCell%WEdge_Area(i+1,j)+1.d-30)/PGrid%dy(i+1,j)
             else
               uright=0.d0
             end if
!             TVar%u(i,j)=0.5d0*(uleft+uright)
!             print*,'Interpolation var'
!             print*,i,j
!             print*,uleft,uright
!             pause'pause to test'
           end if
           ! from momentum exchange cell to normal cut cell
           if(VCellO%MoExCell(i,j)==1.and.VCell%MoExCell(i,j)==0) then
             ii = VCellO%MsCe(i,j,1)
             jj = VCellO%MsCe(i,j,2)
           !  TVar%v(ii,jj) = TVar%v(ii,jj)*(VCell%vof(i,j)+VCell%vof(ii,jj))/  &
           !                  (VCell%vof(ii,jj)+VCell%vof(i,j)*VCell1%delh(i,j)/&
           !                  VCell1%delh(ii,jj))
           !  vfa = vb+(TVar%v(i,j)-vb)*VCell%delh(i,j)/VCell%delh(ii,jj)
           !  TVar%v(ii,jj) = (TVar%v(ii,jj)*VCell%vof(ii,jj)+TVar%v(i,j)*      &
           !                  VCell%vof(i,j)-vfa*VCell%vof(i,j))/VCell%vof(ii,jj)
           !  TVar%v(i,j) = vfa
           end if
           ! from normal cut cell to momentum exchange cell
           if(VCell%MoExCell(i,j)==1.and.VCellO%MoExCell(i,j)==0) then
             ii = VCell%MsCe(i,j,1)
             jj = VCell%MsCe(i,j,2)
           !  TVar%v(ii,jj) = (TVar%v(i,j)*VCell%vof(i,j)+TVar%v(ii,jj)*        &
           !                 VCell%vof(ii,jj))/(VCell%vof(i,j)+VCell%vof(ii,jj))
           !  TVar%v(i,j) = TVar%v(ii,jj)
           end if
           ! for V-velocity (new velocity cell)
           if(VCellO%VofS(i,j)>1.d0-epsi.and.VCell%VofS(i,j)<1.d0-epsi) then
             ii = VCell%MsCe(i,j,1)
             jj = VCell%MsCe(i,j,2)
             TVar%v(i,j) = vb !TVar%v(ii,jj)!
             if(PCell%Nedge_Area(i,j)>epsi) then

               vbot=((-TVar%u(i,j)*PCell%EEdge_Area(i,j)                       &
                      +TVar%u(i-1,j)*PCell%WEdge_Area(i,j))*PGrid%dy(i,j)      &
                      +TVar%v(i,j-1)*PCell%SEdge_Area(i,j)*PGrid%dx(i,j))      &
                      /PGrid%dx(i,j)/(PCell%NEdge_Area(i,j)+1.d-30)
             else
               vbot=vb
             end if
             if(PCell%SEdge_Area(i,j+1)>epsi) then
               vtop=((TVar%u(i,j+1)*PCell%EEdge_Area(i,j+1)                    &
                     -TVar%u(i-1,j+1)*PCell%WEdge_Area(i,j+1))*PGrid%dy(i,j+1) &
                     +TVar%v(i,j+1)*PCell%NEdge_Area(i,j+1)*PGrid%dx(i,j+1))   &
                     /PGrid%dx(i,j+1)/(PCell%SEdge_Area(i,j+1)+1.d-30)
             else
               vbot=vb
             end if
!             TVar%v(i,j)=0.5d0*(vbot+vtop)
           end if
           if(PCell%NEdge_Area(i,j)>epsi.and.PCellO%NEdge_Area(i,j)>epsi) then
             if(PCell%NEdge_Area(i,j)>0.5d0.or.PCellO%NEdge_Area(i,j)>0.5d0) then
               if(PCellO%NEdge_Area(i,j)/PCell%NEdge_Area(i,j)<0.6d0) then
                 TVar%v(i,j)=TVar%v(i,j)*PCellO%NEdge_Area(i,j)/PCell%NEdge_Area(i,j)
!                 print*,'this for interpolation VVVVVVVVVVVVVVV'
!                 print*,i,j
!                 print*,PCellO%NEdge_Area(i,j),PCell%NEdge_Area(i,j)
!                 print*,PCellO%NEdge_Area(i,j)/PCell%NEdge_Area(i,j)
!                 print*,TVar%v(i,j)
!                 print*,'*************************'
               end if
             end if
           end if
           if(PCell%EEdge_Area(i,j)>epsi.and.PCellO%EEdge_Area(i,j)>epsi) then
             if(PCell%EEdge_Area(i,j)>0.5d0.or.PCellO%EEdge_Area(i,j)>0.5d0) then
               if(PCellO%EEdge_Area(i,j)/PCell%EEdge_Area(i,j)<0.6) then
                 TVar%u(i,j)=TVar%u(i,j)*PCellO%EEdge_Area(i,j)/PCell%EEdge_Area(i,j)
!                 print*,'this for interpolation UUUUUUUUUUUUUUU'
!                 print*,i,j
!                 print*,PCellO%EEdge_Area(i,j),PCell%EEdge_Area(i,j)
!                 print*,PCellO%EEdge_Area(i,j)/PCell%EEdge_Area(i,j)
!                 print*,TVar%u(i,j)
!                 print*,'*************************'
               end if
             end if
           end if
         end do
       end do
    END SUBROUTINE InterNewVar

    SUBROUTINE ReadFileConvergence(filename,Time,TNorm,BoomCase)
      IMPLICIT NONE
      CHARACTER(LEN=80),INTENT(IN)         :: filename
      TYPE(SolverTime),INTENT(INOUT)       :: Time
      TYPE(SolverConvergence),INTENT(INOUT):: TNorm
      TYPE(SolidObject),INTENT(INOUT)      :: BoomCase
      INTEGER(KIND=it4b)                   :: i,j
      INTEGER(it8b) :: IttRun
      call getSolverVariables(IttRun_=IttRun)
      open(unit=5,file=filename,action='read')
      do i=1,IttRun
        read(5,76) j,Time%NondiT,TNorm%N1,TNorm%N2,TNorm%Ninf,         &
                     TNorm%N1c, TNorm%N2c,TNorm%Ninfc,BoomCase%Posp%y, &
                     BoomCase%PospO%y,BoomCase%vs,BoomCase%asy
      end do
      close(5)
  76  format(I10,11(f24.14))
    END SUBROUTINE ReadFileConvergence


    SUBROUTINE PrintWaterHeight(NonTime,I1,I2,I3,I4,PGrid,PCell)
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(IN):: NonTime
      INTEGER(kind=it4b),INTENT(IN):: I1,I2,I3,I4
      TYPE(Grid),INTENT(IN):: PGrid
      TYPE(Cell),INTENT(IN):: PCell
      INTEGER(kind=it4b):: j
      REAL(KIND=dp):: H1,H2,H3,H4
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      CHARACTER*70 :: dir
      call getDir(dir)
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      do j = 1,Jsize
        if(PCell%vof(I1,j)>epsi.and.PCell%vof(I1,j)<1.d0-epsi) then
          H1 = PGrid%y(I1,j)+(PCell%vof(I1,j)-0.5d0)*PGrid%dy(I1,j)
        end if
        if(PCell%vof(I2,j)>epsi.and.PCell%vof(I2,j)<1.d0-epsi) then
          H2 = PGrid%y(I2,j)+(PCell%vof(I2,j)-0.5d0)*PGrid%dy(I2,j)
        end if
        if(PCell%vof(I3,j)>epsi.and.PCell%vof(I3,j)<1.d0-epsi) then
          H3 = PGrid%y(I3,j)+(PCell%vof(I3,j)-0.5d0)*PGrid%dy(I3,j)
        end if
        if(PCell%vof(I4,j)>epsi.and.PCell%vof(I4,j)<1.d0-epsi) then
          H4 = PGrid%y(I4,j)+(PCell%vof(I4,j)-0.5d0)*PGrid%dy(I4,j)
        end if
      end do
      open(unit=5,file=trim(adjustl(dir))//'WaterHeight.dat',access='append')
      write(5,78) NonTime,H1,H2,H3,H4
      close(5)
 78	  format(f15.10,f15.10,f15.10,f15.10,f15.10)
    END SUBROUTINE PrintWaterHeight

    SUBROUTINE PrintWaterWave(TimePrint,PGrid,PCell)
      REAL(KIND=dp),INTENT(IN):: TimePrint
      TYPE(Grid),INTENT(IN):: PGrid
      TYPE(Cell),INTENT(IN):: PCell
      INTEGER(kind=it4b):: i,j
      character(len=250):: curd
      REAL(KIND=dp):: VolPrint,XPrint,HPrint,tol,Hmax
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      CHARACTER*70 :: dir
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      call getDir(dir)

      tol=1.d-24
      write(curd,'(f15.8)') TimePrint
      open(unit=5,file=trim(adjustl(dir))//'WaveProfile_'//trim(curd)//'.dat', &
                                                                action='write')
      if(WavePrint==3) then
        Hmax=2.d0
      else
        Hmax=0.d0
      end if
      do i=1,Isize
        VolPrint=1.d0
        do j=1,Jsize
          if(PCell%vof(i,j)>epsi.and.PCell%vof(i,j)<1.d0-epsi) then
          ! Check for the cell which has liquid volume fraction closing to 0.5
            if(dabs(PCell%vof(i,j)-0.5d0)<VolPrint) then
              VolPrint=dabs(PCell%vof(i,j)-0.5d0)
              XPrint=PGrid%x(i,j)
          ! Calculate the Liquid height
              HPrint=PGrid%y(i,j)-PCell%phi(i,j)/(PCell%ny(i,j)+tol)-          &
                                                              0.05d0/PGrid%Lref
            end if
          end if
        end do
        write(5,100) (XPrint)*100.d0,HPrint*100.d0
      end do
      close(5)
 100  format(f15.10,f15.10)
    END SUBROUTINE PrintWaterWave

    SUBROUTINE PrintWaterPressure(NonTime,J1,J2,J3,J4,Var)
      REAL(KIND=dp),INTENT(IN):: NonTime
      INTEGER(kind=it4b),INTENT(IN):: J1,J2,J3,J4
      TYPE(TVariables),INTENT(IN):: Var
      REAL(KIND=dp):: P1,P2,P3,P4
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      CHARACTER*70 :: dir
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      call getDir(dir)
      P1=Var%p(Isize,J1)
      P2=Var%p(Isize,J2)
      P3=Var%p(Isize,J3)
      p4=Var%p(Isize,J4)
      open(unit=5,file=trim(adjustl(dir))//'WaterPressure.dat',access='append')
      write(5,78) NonTime,P1,P2,P3,P4
      close(5)
 78	  format(f15.10,f15.10,f15.10,f15.10,f15.10)
    END SUBROUTINE PrintWaterPressure

    SUBROUTINE CopyNewCell(TCellCo,TCellTa)
      IMPLICIT NONE
      TYPE(Cell),INTENT(IN):: TCellTa
      TYPE(Cell),INTENT(INOUT):: TCellCo
      INTEGER(kind=it4b):: i,j
      INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize
      call getMeshSizes(ibeg, jbeg, Isize, Jsize)
      do i = ibeg,ibeg+Isize-1
        do j = jbeg,jbeg+Jsize-1
          TCellCo%vof(i,j)=TCellTa%vof(i,j)
          TCellCo%phi(i,j)=TCellTa%phi(i,j)
          TCellCo%nx(i,j)=TCellTa%nx(i,j)
          TCellCo%ny(i,j)=TCellTa%ny(i,j)
          TCellCo%vofS(i,j)=TCellTa%vofS(i,j)
          TCellCo%phiS(i,j)=TCellTa%phiS(i,j)
          TCellCo%nxS(i,j)=TCellTa%nxS(i,j)
          TCellCo%nyS(i,j)=TCellTa%nyS(i,j)

          TCellCo%EEdge_Area(i,j)=TCellTa%EEdge_Area(i,j)
          TCellCo%WEdge_Area(i,j)=TCellTa%WEdge_Area(i,j)
          TCellCo%NEdge_Area(i,j)=TCellTa%NEdge_Area(i,j)
          TCellCo%SEdge_Area(i,j)=TCellTa%SEdge_Area(i,j)
          if(allocated(TCellCo%MoExCell)) then
            TCellCo%MoExCell(i,j)=TCellTa%MoExCell(i,j)
          end if
          if(allocated(TCellCo%Cell_Cent)) then
            TCellCo%Cell_Cent(i,j,:)=TCellTa%Cell_Cent(i,j,:)
          end if
          if(allocated(TCellCo%EtaE))TCellCo%EtaE(i,j)=TCellTa%EtaE(i,j)
          if(allocated(TCellCo%EtaE))TCellCo%EtaN(i,j)=TCellTa%EtaN(i,j)
          if(allocated(TCellCo%AlE))TCellCo%AlE(i,j)=TCellTa%AlE(i,j)
          if(allocated(TCellCo%AlN))TCellCo%AlN(i,j)=TCellTa%AlN(i,j)
          if(allocated(TCellCo%DAlE))TCellCo%DAlE(i,j)=TCellTa%DAlE(i,j)
          if(allocated(TCellCo%DAlN))TCellCo%DAlN(i,j)=TCellTa%DAlN(i,j)
          if(allocated(TCellCo%SxE))TCellCo%SxE(i,j)=TCellTa%SxE(i,j)
          if(allocated(TCellCo%SyN))TCellCo%SyN(i,j)=TCellTa%SyN(i,j)
          if(allocated(TCellCo%FCE))TCellCo%FCE(i,j,:)=TCellTa%FCE(i,j,:)
          if(allocated(TCellCo%FCN))TCellCo%FCN(i,j,:)=TCellTa%FCN(i,j,:)
          if(allocated(TCellCo%WlLh))TCellCo%WlLh(i,j)=TCellTa%WlLh(i,j)
          if(allocated(TCellCo%delh))TCellCo%delh(i,j)=TCellTa%delh(i,j)
          if(allocated(TCellCo%MsCe)) TCellCo%MsCe(i,j,:)=TCellTa%MsCe(i,j,:)
          TCellCo%ExtCell=TCellTa%ExtCell
        end do
      end do
      do i = ibeg-1,ibeg+Isize-1+1
        do j = jbeg-1,jbeg+Jsize-1+1
          TCellCo%EEdge_Area(i,j)=TCellTa%EEdge_Area(i,j)
          TCellCo%WEdge_Area(i,j)=TCellTa%WEdge_Area(i,j)
          TCellCo%NEdge_Area(i,j)=TCellTa%NEdge_Area(i,j)
          TCellCo%SEdge_Area(i,j)=TCellTa%SEdge_Area(i,j)
        end do
      end do
    END SUBROUTINE CopyNewCell
!    SUBROUTINE WaveBoundaryCondition(PGrid,Vari,Time)
!      TYPE(Grid),INTENT(IN):: PGrid
!      TYPE(TVariables),INTENT(INOUT):: Vari
!      TYPE(SolverTime),INTENT(IN):: Time
!      INTEGER(kind=it8b):: i,j
!      REAL(KIND=dp):: etau,etav
!      etau=amp0*dsin(2.d0*pi/Lamdaw*(-cw0*Time%NondiT))
!      etav=amp0*dsin(2.d0*pi/Lamdaw*(-PGrid%dx(1,1)/2.d0-cw0*Time%NondiT))
!      do j=jbeg,Jsize+jbeg-1
!      ! Left inlet wall
!        Vari%p(ibeg-ight,j)=Vari%p(ibeg,j) !Vari%Pint/(Vari%Pref)
!        Vari%t(ibeg-ight,j)=Vari%t(ibeg,j)
!        if(PGrid%y(1,j)-PGrid%dy(1,j)/2.d0<Hw+etau) then
!        ! for water entry velocity
!          Vari%u(ibeg-ight,j)=2.d0*pi/Lamdaw*cw0*etau*                         &
!                 dcosh(2.d0*pi*PGrid%y(1,j))/dsinh(2.d0*pi*Hw)
!        else
!        ! for gas entry velocity
!          Vari%u(ibeg-ight,j)=-2.d0*pi/Lamdaw*cw0*etau*                        &
!                 dcosh(2.d0*pi*(HChannel-PGrid%y(1,j)))/dsinh(2.d0*pi*Ha)
!        end if
!        if(PGrid%y(1,j)-PGrid%dy(1,j)/2.d0<Hw+etav) then
!          Vari%v(ibeg-ight,j)=-Amp0*2.d0*pi/cw0*                               &
!                dcos(2.d0*pi/Lamdaw*(-PGrid%dx(1,1)/2.d0-cw0*Time%NondiT))*    &
!                dsinh(2.d0*pi*PGrid%y(1,j))/dsinh(2.d0*pi*Hw)
!        else
!          Vari%v(ibeg-ight,j)=Amp0*2.d0*pi/cw0*                                &
!                dcos(2.d0*pi/Lamdaw*(-PGrid%dx(1,1)/2.d0-cw0*Time%NondiT))*    &
!                dsinh(2.d0*pi*(HChannel-PGrid%y(1,j)))/dsinh(2.d0*pi*Ha)
!        end if
!      ! Right Wall non-Slip wall
!        Vari%p(Isize+ibeg+ight-1,j)=Vari%p(Isize+ibeg-1,j)
!        Vari%t(Isize+ibeg+ight-1,j)=Vari%t(ibeg+Isize-1,j)
!        Vari%u(Isize+ibeg+ight-1,j)=0.d0!-Vari%u(Isize+ibeg-1,j) !Vari%Uint/Vari%Uref
!        Vari%u(Isize+ibeg-1,j)=0.d0
!        Vari%v(Isize+ibeg+ight-1,j)=0.d0-Vari%v(Isize+ibeg-1,j) !Vari%Vint/Vari%Uref
!      end do
!      do i = ibeg,Isize+ibeg-1
!      ! Bottom wall non-slip wall
!        Vari%p(i,jbeg-jght)=Vari%p(i,jbeg)
!        Vari%t(i,jbeg-jght)=Vari%t(i,jbeg)
!        Vari%u(i,jbeg-jght)=Vari%u(i,jbeg)
!        Vari%v(i,jbeg-jght)=0.d0 !Vari%Vint/Vari%Uref
!      ! Open air
!        Vari%p(i,Jsize+jbeg+jght-1)=0.d0
!        Vari%t(i,jbeg+Jsize+jght-1)=Vari%t(i,jbeg+Jsize-1)
!        Vari%u(i,Jsize+jbeg+jght-1)=Vari%u(i,Jsize+jbeg-1)
!     !   Vari%v(i,Jsize+jbeg-1) = Vari%v(i,Jsize+jbeg-2)
!     !   Vari%v(i,Jsize+jbeg+jght-1) = Vari%v(i,Jsize+jbeg-1)
!      end do
!    END SUBROUTINE WaveBoundaryCondition
END MODULE Solver
