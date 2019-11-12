    !*******************************************************
    !
    !              Open air
    !     |                        |
    !     |Air                     |
    !     |_____                   |
    !     |     |                  |
    ! Wall|Water|                  |Wall
    !     |     |                  |
    !     |     |                  |
    !     |_____|__________________|
    !                Wall
    !*******************************************************
Program Main
    USE PrecisionVar
    USE Mesh, ONLY : TsimcoMesh, TPoint
    USE Cutcell, ONLY : Grid_Preprocess, NewCellFace
    USE Clsvof, ONLY : SolidObject, Initial_Clsvof
    USE StateVariables, ONLY : TVariables, TWave, setSolverVariables
    USE Constants, ONLY : g, pi, nuw, row
    USE PrintResult, ONLY : Print_Result_Tecplot_PCent, Print_Result_Tecplot_UCent, Print_Result_Tecplot_VCent, setDir
    USE MPI, ONLY : MPI_Initial
    USE Solver, ONLY : IterationSolution, IterationSolution2
    USE Particles, ONLY : TParticle
    USE BoundaryFunction, ONLY : BCBase
    USE BoundaryFunction2!, ONLY : BCBase2, BCUW

    IMPLICIT NONE
    TYPE(TPoint):: ReS,ReE, Start_Point, End_Point
    TYPE(TVariables)  :: Var
    TYPE(TWave)       :: wave
    TYPE(TParticle)   :: TraPar
    TYPE(SolidObject) :: BoomCase
    TYPE(BCBase)      :: BCp, BCu, BCv, BCVof, BCLvs
    TYPE(BCBase2)     :: BCp2, BCu2, BCv2, BCVof2, BCLvs2
    INTEGER(kind=it4b):: Irec,Jrec,NI,NJ,iprint
    REAL(dp) :: zp,UParInlet,HParInlet,DParInlet,rop,gx,gy, Lref, Roref, Rey
    REAL(KIND=dp):: vel,Uref,Vint, ha, hw, UwInlet, UgInlet, nuref
    REAL(dp) :: t0, cw0, Amp0, Depthw, Lamdaw, twp, HChannel, LChannel, kw, omew, HDomain
    REAL(dp), DIMENSION(:), ALLOCATABLE :: Constin
    TYPE(TsimcoMesh) :: simcomesh
    INTEGER(it4b) :: Isize, Jsize, NParInlet, IParInlet
    LOGICAL :: RunAgain, ICorProb
    INTEGER(it8b) :: IttRun 
    CHARACTER*70 :: dir
    open(unit=5,file='input.dat',action='read')
    read(5,*)
    read(5,*) Isize, Jsize, Irec, Jrec, Rey, Hw, iprint
    read(5,*)
    read(5,'(a70)') dir
    close(5)
    call setDir(dir)
    Isize=500
    Jsize=160
    NI=Isize+1
    NJ=Jsize+1
    call MPI_Initial
    simcomesh = TsimcoMesh(Isize, Jsize)
    Var       = Tvariables()

    !Particles
    IParInlet=50
    NParInlet=5
    TraPar = TParticle(0, NParInlet, IParInlet) ! Zero particles

    Lref=1.d0 !Hw
    Uref=1.d0
    Roref=row
    nuref=nuw
    rop=711.d0
    Hw=6.0d0/Lref !DWedge
    HDomain=8.d0
    HChannel=80.d0/Lref
    LChannel=25.d0/Lref
    Ha=HDomain-Hw
    Depthw=HChannel-Ha
  ! for sinusoidal wave
    twp=1.131371d0
    Lamdaw=g*twp**2.d0/2.d0/pi
    cw0=dsqrt(g*Lamdaw/2.d0/pi)
    Lamdaw=Lamdaw/Lref
    cw0=cw0/Uref
    twp=twp/(Lref/Uref)
    kw=2.d0*pi/Lamdaw
    Amp0=0.1d0/Lref/2.d0
    omew=2.d0*pi/twp
    Start_Point = TPoint(0.d0/Lref, (HChannel-HDomain)/Lref) !-Hw/Lref
    End_Point = TPoint(LChannel/Lref, HChannel/Lref)
    BoomCase%Posp%x=15.d0/Lref
    BoomCase%Posp%y=Depthw
    BoomCase%Dobj=0.8d0/Lref
    BoomCase%Wobj=0.16d0/Lref
    BoomCase%XBar1=BoomCase%Posp%x-BoomCase%Wobj/2.d0
    BoomCase%XBar2=BoomCase%Posp%x+BoomCase%Wobj/2.d0
    BoomCase%LBar=1.5d0/Lref
    BoomCase%YBar=BoomCase%Posp%y-dsqrt((BoomCase%Dobj/2.d0)**2.d0-                &
                 (BoomCase%Wobj/2.d0)**2.d0)-BoomCase%LBar
    BoomCase%Mobj=(pi/4.d0*(BoomCase%Dobj)**2.d0+BoomCase%Wobj*BoomCase%LBar)*0.5d0*rop/Roref

    wave = TWave(t0, cw0, Amp0, Depthw, Lamdaw, twp, HChannel, LChannel, kw, omew, HDomain)
    UwInlet=1.d0/Uref
    UgInlet=3.d0/Uref
    UParInlet=4.d0/Uref
    HParInlet=1.5d0/Lref
    DParInlet=1.d-2/Lref
    ReS = TPoint(13.5d0/Lref, 0.2d0/Lref)
    ReE = TPoint(16.0d0/Lref, 0.8d0/Lref)
    zp=1.d-2/Lref
    Irec=375
    Jrec=60
    Vint=0.d0
    vel=0.d0
    gx=0.d0
    gy=g
    RunAgain=.FALSE.
!   define corner problem
    ICorProb=.TRUE.
    IttRun=600
    CALL setSolverVariables(RunAgain, ICorProb, IttRun)
!    call Initial_Grid(Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,PGrid,Lref,0)
!    call InitialUVGrid(PGrid,UGrid,0,Lref)
!    call InitialUVGrid(PGrid,VGrid,1,Lref)
    ! Using TsimcoMesh
    call simcomesh%Initial_Grid2(Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,Lref,0)
    call simcomesh%InitialUVGrid2(Lref, Lref)
    call simcomesh%HYPRE_CreateGrid2()
    ! Contruct the layer of boundary points
    call BCp%Construct(Isize,Jsize)
    BCp2 = BCBase2(Isize, Jsize)
    call BCu%Construct(Isize,Jsize)
    BCu2 = BCBase2(Isize, Jsize)
    call BCv%Construct(Isize,Jsize)
    BCv2 = BCBase2(Isize, Jsize)
    call BCVof%Construct(Isize,Jsize)
    BCVof2 = BCBase2(Isize, Jsize)
    call BCLvs%Construct(Isize,Jsize) 
    BCLvs2 = BCBase2(Isize, Jsize)
    ! Set up the flag for Dirichlet or Neumann boundary type.
    ! The Boundary condition is decided by the user. 
    call BCp%SetDN(1,0,1,1) 
    call BCu%SetDN(0,1,0,0)
    call BCv%SetDN(0,1,0,0)
    call BCvof%SetDN(1,1,0,0)
    call BCp2%SetDN(1,0,1,1) 
    call BCu2%SetDN(0,1,0,0)
    call BCv2%SetDN(0,1,0,0)
    call BCvof2%SetDN(1,1,0,0)
    ! Note that the Neumann boundary condition will reduce the computational time 
    ! for computing level set function at boundary
    call BClvs%SetDN(1,1,1,1)
    call BClvs2%SetDN(1,1,1,1)
    ! Set constant for computing value at boundary
    ! pressure
    allocate(Constin(4))
    Constin(:)=0.d0
    call BCp%SetConstant(Constin)
    call BCp2%SetConstant(Constin)
    deallocate(Constin)
    allocate(Constin(11))
    Constin(1:4) = 0.d0
    ! Set constant for wave 
    Constin(5) = Amp0
    Constin(6) = Kw
    Constin(7) = cw0
    Constin(8) = UwInlet
    Constin(9) = UGInlet
    Constin(10) = HChannel
    Constin(11) = Depthw
    call BCu%SetConstant(Constin)
    call BCv%SetConstant(Constin)
    call BCVof%SetConstant(Constin)
    call BCLvS%SetConstant(Constin)
    call BCu2%SetConstant(Constin)
    call BCv2%SetConstant(Constin)
    call BCVof2%SetConstant(Constin)
    call BCLvS2%SetConstant(Constin)
    BCu2%north=>BCuN
    BCu2%south=>BCuS
    BCu2%east =>BCuE 
    BCu2%west =>BCuW
    BCv2%north=>BCvN
    BCv2%south=>BCvS
    BCv2%east =>BCvE 
    BCv2%west =>BCvW
    BCp2%north=>BCpN
    BCp2%south=>BCpS
    BCp2%east =>BCpE 
    BCp2%west =>BCpW
    deallocate(Constin)        
!    call HYPRE_CreateGrid(PGrid)
    call Initial_Clsvof(simcomesh%PGrid,simcomesh%PCell,BoomCase, wave)
    call Initial_Clsvof(simcomesh%UGrid,simcomesh%UCell,BoomCase, wave)
    call Initial_Clsvof(simcomesh%VGrid,simcomesh%VCell,BoomCase, wave)
    call Grid_Preprocess(simcomesh,Var,int8(1))
    call Var%Initialize(wave, simcomesh, simcomesh%PCell,simcomesh%PGrid,vel,Vint,0.d0,300.d0,Uref,300.d0,&
         &              Roref,Lref, nuref, ha, UwInlet, UgInlet)
    call TraPar%InitializeParticles(simcomesh%PGrid,Var, UParInlet, HParInlet, DParInlet, zp, rop, gx, gy)
  !  call Print_Result_Tecplot_PCent(PGrid,Var,PCell,TraPar,INT8(0),1)
    call Print_Result_Tecplot_UCent(simcomesh%UGrid,Var,simcomesh%UCell,INT8(0))
    call Print_Result_Tecplot_VCent(simcomesh%VGrid,Var,simcomesh%VCell,INT8(0))
    call NewCellFace(simcomesh)
    call IterationSolution2(simcomesh, Var, wave, TraPar, BoomCase, BCp2, BCu2, BCv2, BCvof2, BCLvs2, 50)
end program main
