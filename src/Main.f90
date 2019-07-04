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
    USE Mesh, ONLY : TsimcoMesh, getMeshSizes, Point
    USE Cutcell, ONLY : Grid_Preprocess, NewCellFace
    USE Clsvof, ONLY : SolidObject, Initial_Clsvof
    USE StateVariables!, ONLY : TVariables, TWave
    USE Constants, ONLY : g, pi, nuw, row
    USE PrintResult, ONLY : Print_Result_Tecplot_PCent, Print_Result_Tecplot_UCent, Print_Result_Tecplot_VCent
    USE MPI, ONLY : MPI_Initial
    USE Solver, ONLY : IterationSolution
    USE Particles, ONLY : Particle, InitializeParticles
    IMPLICIT NONE
    TYPE(Point):: ReS,ReE, Start_Point, End_Point
    TYPE(TVariables):: Var
    TYPE(TWave)     :: wave
    TYPE(Particle):: TraPar
    TYPE(SolidObject):: BoomCase
    INTEGER(kind=it4b):: Irec,Jrec,NI,NJ,iprint
    REAL(KIND=dp):: vel,Uref,Vint
    REAL(dp) :: t0, cw0, Amp0, Depthw, Lamdaw, twp, HChannel, LChannel, kw, omew, HDomain
    TYPE(TsimcoMesh) :: simcomesh
    INTEGER(it4b) :: ibeg, jbeg, Isize, Jsize, ight, jght
    call getMeshSizes(ibeg, jbeg, ighte=ight, jghte=jght)
    open(unit=5,file='input.dat',action='read')
    read(5,*),
    read(5,*), Isize, Jsize, Irec, Jrec, Rey, Hw, iprint
    read(5,*),
    read(5,'(a70)'), dir
    close(5)
    Isize=500
    Jsize=160
    NI=Isize+1
    NJ=Jsize+1
    call MPI_Initial
    simcomesh = TsimcoMesh(Isize, Jsize)
    Var       = Tvariables()

    allocate(TraPar%dp(10000))
    allocate(TraPar%Posp(10000))
    allocate(TraPar%uvp(10000))
    TraPar%np=0
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
 !   Start_Point%x=-6.5/Lref
    Start_Point%x=0.d0/Lref
    Start_point%y=(HChannel-HDomain)/Lref !-Hw/Lref
 !   End_Point%x=6.5d0/Lref
    End_Point%x=LChannel/Lref
    End_Point%y=HChannel/Lref
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
    IParInlet=50
    NParInlet=5
    ReS%x=13.5d0/Lref
    ReS%y=0.2d0/Lref
    ReE%x=16.0d0/Lref
    ReE%y=0.8d0/Lref
    zp=1.d-2/Lref
    Irec=375
    Jrec=60
    Vint=0.d0
    vel=0.d0
    gx=0.d0
    gy=g
    VofInlet=1.d0
    RunAgain=.FALSE.
!   define corner problem
    ICorProb=.TRUE.
    IttRun=600
!    call Initial_Grid(Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,PGrid,Lref,0)
!    call InitialUVGrid(PGrid,UGrid,0,Lref)
!    call InitialUVGrid(PGrid,VGrid,1,Lref)
    ! Using TsimcoMesh
    call simcomesh%InitialUVGrid2(Lref, Lref)
    call simcomesh%Initial_Grid2(Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,Lref,0)
    call simcomesh%HYPRE_CreateGrid2()
    !
!    call HYPRE_CreateGrid(PGrid)
    call Initial_Clsvof(simcomesh%PGrid,simcomesh%PCell,BoomCase, wave)
    call Initial_Clsvof(simcomesh%UGrid,simcomesh%UCell,BoomCase, wave)
    call Initial_Clsvof(simcomesh%VGrid,simcomesh%VCell,BoomCase, wave)
    call Grid_Preprocess(simcomesh,Var,int8(1))
    call Var%Initialize(wave, simcomesh, simcomesh%PCell,simcomesh%PGrid,vel,Vint,0.d0,300.d0,Uref,300.d0,Roref,Lref)
    call InitializeParticles(simcomesh%PGrid,Var,TraPar)
  !  call Print_Result_Tecplot_PCent(PGrid,Var,PCell,TraPar,INT8(0),1)
    call Print_Result_Tecplot_UCent(simcomesh%UGrid,Var,simcomesh%UCell,INT8(0))
    call Print_Result_Tecplot_VCent(simcomesh%VGrid,Var,simcomesh%VCell,INT8(0))
    call NewCellFace(simcomesh)
    call IterationSolution(simcomesh, Var,wave, TraPar,BoomCase,50)
    pause
end program main
