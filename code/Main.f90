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
    USE Mesh
    USE Cutcell
    USE Clsvof
    USE StateVariables
    USE PrintResult
    USE MPI
    USE Solver
    USE Particles
    IMPLICIT NONE
    TYPE(Grid):: UGrid,VGrid,PGrid
    TYPE(Cell):: UCell,VCell,PCell
    TYPE(Point):: ReS,ReE
    TYPE(Variables):: Var
    TYPE(Particle):: TraPar
    TYPE(SolidObject):: BoomCase
    INTEGER(kind=it4b):: Irec,Jrec,NI,NJ,iprint
    REAL(KIND=dp):: vel,Uref,Vint
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
    allocate(UGrid%x(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UGrid%y(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VGrid%x(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VGrid%y(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(PGrid%x(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(PGrid%y(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UGrid%dx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UGrid%dy(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VGrid%dx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VGrid%dy(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(PGrid%dx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(PGrid%dy(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))

    allocate(UCell%vof(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%phi(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%nx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%ny(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%vofS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%phiS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%nxS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%nyS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%EEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
    allocate(UCell%WEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
    allocate(UCell%NEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
    allocate(UCell%SEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
    allocate(UCell%MoExCell(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(UCell%EtaE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%EtaN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%DAlE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%DAlN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%AlE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%AlN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%SxE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%SyN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%FCE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
    allocate(UCell%FCN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
    allocate(UCell%WlLh(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(UCell%delh(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(UCell%PosNu(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(UCell%Cell_Cent(Isize,Jsize,2))
    allocate(UCell%MsCe(Isize,Jsize,2))

    allocate(VCell%vof(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%phi(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%nx(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%ny(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%vofS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%phiS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%nxS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%nyS(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%EEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
    allocate(VCell%WEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
    allocate(VCell%NEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
    allocate(VCell%SEdge_Area(ibeg-1:Isize+ibeg-1+1,jbeg-1:Jsize+jbeg-1+1))
    allocate(VCell%MoExCell(ibeg:Isize+ibeg-1,jbeg:Jsize+jbeg-1))
    allocate(VCell%EtaE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%EtaN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%DAlE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%DAlN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%AlE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%AlN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%SxE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%SyN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%FCE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
    allocate(VCell%FCN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
    allocate(VCell%WlLh(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(VCell%delh(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(VCell%PosNu(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(VCell%Cell_Cent(Isize,Jsize,2))
    allocate(VCell%MsCe(Isize,Jsize,2))

    allocate(PCell%vof(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%phi(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%nx(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%ny(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%vofS(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%phiS(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%nxS(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%nyS(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%EEdge_Area(ibeg-1:Isize+1,jbeg-1:Jsize+1))
    allocate(PCell%WEdge_Area(ibeg-1:Isize+1,jbeg-1:Jsize+1))
    allocate(PCell%NEdge_Area(ibeg-1:Isize+1,jbeg-1:Jsize+1))
    allocate(PCell%SEdge_Area(ibeg-1:Isize+1,jbeg-1:Jsize+1))
    allocate(PCell%FCE(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
    allocate(PCell%FCN(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght,2))
    allocate(PCell%WlLh(ibeg:Isize,jbeg:Jsize))
    allocate(PCell%PosNu(ibeg-ight:ibeg+Isize-1+ight,jbeg-jght:jbeg+Jsize-1+jght))
    allocate(PCell%Cell_Cent(Isize,Jsize,2))

    allocate(Var%u(ibeg-ight:Isize-ibeg+ight+1,jbeg-jght:Jsize-jbeg+jght+1))
    allocate(Var%v(ibeg-ight:Isize-ibeg+ight+1,jbeg-jght:Jsize-jbeg+jght+1))
    allocate(Var%p(ibeg-ight:Isize-ibeg+ight+1,jbeg-jght:Jsize-jbeg+jght+1))
    allocate(Var%Gpu(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(Var%Gpv(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(Var%ures(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(Var%vres(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(Var%pres(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(Var%mres(ibeg:ibeg+Isize-1,jbeg:jbeg+Jsize-1))
    allocate(Var%t(ibeg-ight:Isize-ibeg+ight+1,jbeg-jght:Jsize-jbeg+jght+1))
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
    call Initial_Grid(Start_Point,End_Point,ReS,ReE,NI,NJ,Irec,Jrec,PGrid,Lref,0)
    call InitialUVGrid(PGrid,UGrid,0,Lref)
    call InitialUVGrid(PGrid,VGrid,1,Lref)
    call MPI_Initial
    call HYPRE_CreateGrid(PGrid)
    call Initial_Clsvof(PGrid,PCell,BoomCase)
    call Initial_Clsvof(UGrid,UCell,BoomCase)
    call Initial_Clsvof(VGrid,VCell,BoomCase)
    call Grid_Preprocess(PGrid,UGrid,VGrid,PCell,UCell,VCell,Var,int8(1))
    call Initial_Var(PCell,PGrid,Var,vel,Vint,0.d0,300.d0,Uref,300.d0,Roref,Lref)
    call InitializeParticles(PGrid,Var,TraPar)
  !  call Print_Result_Tecplot_PCent(PGrid,Var,PCell,TraPar,INT8(0),1)
    call Print_Result_Tecplot_UCent(UGrid,Var,UCell,INT8(0))
    call Print_Result_Tecplot_VCent(VGrid,Var,VCell,INT8(0))
    call NewCellFace(PCell,UCell,VCell,PGrid,UGrid,VGrid)
    call IterationSolution(PGrid,UGrid,VGrid,PCell,UCell,VCell,Var,TraPar,BoomCase,50)
    pause
end program main
