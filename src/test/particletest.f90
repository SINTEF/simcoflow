!FIXME: This has too many dependencies for now
Program particletest
   USE PrecisionVar
   USE Particles, ONLY : TParticle, ParticlePosition
   USE StateVariables, ONLY :TVariables
   USE Mesh, ONLY : TsimcoMesh, TPoint
   USE Constants, ONLY : g
   IMPLICIT NONE
   TYPE(TParticle)  :: TraPar
   TYPE(TVariables) :: Var
   TYPE(TPoint)     :: point, point2, ReS, ReE
   TYPE(TsimcoMesh) :: mesh
   INTEGER(it4b) :: ii, jj, IParInlet, NParInlet
   !Need to create a mesh first to set Isize and Jsize used all over
   mesh = TsimcoMesh(2,3)
   !Particles
   IParInlet=50
   NParInlet=5
   TraPar = TParticle(0, NParInlet, IParInlet) ! Zero particles
   Var    = Tvariables()
   call TraPar%InitializeParticles(mesh%PGrid,Var, 1.d0, 1.d0, 1.d-2, 1.d-2, 711.d0, 0d0, g)
   point = TPoint(1.01d0,1.01d0)
   point2 = TPoint(0.0d0, 0.0d0)
   ReS = TPoint(13.5d0, 0.2d0)
   ReE = TPoint(16.0d0, 0.8d0)
   call mesh%Initial_Grid2(point2,point,ReS,ReE,3,4,1,2,1.d0,0)
   call ParticlePosition(point, mesh%PGrid, ii, jj)
   IF (.NOT. ii==2 .AND. .NOT. jj == 3) CALL exit(-1)
   call ParticlePosition(point2, mesh%PGrid, ii, jj)
   IF (.NOT. ii==1 .AND. .NOT. jj == 1) CALL exit(-1)

END PROGRAM particletest
