!FIXME: This has too many dependencies for now
Program particletest
   USE PrecisionVar
   USE Particles, ONLY : TParticle, ParticlePosition
   USE StateVariables, ONLY :TVariables
   USE Mesh, ONLY : TsimcoMesh, TPoint
   USE Constants, ONLY : g,row,roa
   IMPLICIT NONE
   TYPE(TParticle)  :: TraPar
   TYPE(TVariables) :: Var
   TYPE(TPoint)     :: point, point2, pointP, ReS, ReE
   TYPE(TsimcoMesh) :: mesh
   INTEGER(it4b) :: i, ii, jj, IParInlet, NParInlet
   REAL(dp)      :: dt, gx, gy,time,romix,gmix
   ! Need to create a mesh first to set Isize and Jsize used all over
   mesh = TsimcoMesh(10,12)
   ! Particles
   IParInlet=50
   NParInlet=5
   Var = Tvariables()
   Var%p(:,:)=0.d0
   Var%u(:,:)=0.d0
   Var%v(:,:)=0.d0
   mesh%PCell%vof(:,:)=1.d0
   mesh%PCell%vofS(:,:)=0.d0
   mesh%PCell%nx(:,:)=1.d0
   mesh%PCell%ny(:,:)=0.d0
   mesh%PCell%nxS(:,:)=1.d0
   mesh%PCell%nyS(:,:)=0.d0
   
   point = TPoint(0.d0,0.d0)
   point2 = TPoint(10.d0,10.d0)
   pointP = TPoint(4.d0,4.d0)
   ReS = TPoint(2.d0,2.d0)
   ReE = TPoint(6.d0,6.d0)
   gx=0.d0
   gy=g
   
   call mesh%Initial_Grid2(point,point2,ReS,ReE,11,13,1,2,1.d0,0)
   ! Intialize particle
   TraPar%np=1
   TraPar%zp  = 1.d0
   TraPar%rop = 500.d0
   TraPar%gx  = gx
   TraPar%gy  = gy
   TraPar%dp(TraPar%np)=1.d-3
   TraPar%Posp(TraPar%np)%x=pointp%x
   TraPar%Posp(Trapar%np)%y=pointp%y
   TraPar%uvp(TraPar%np)%u=0.d0
   TraPar%uvp(TraPar%np)%v=0.d0
   time = 0.d0
   do i=1,100
     dt = 0.5d0*dmin1(mesh%PGrid%dx(1,1)/(TraPar%uvp(TraPar%np)%u+1.d-30+    &
                      dsqrt(mesh%PGrid%dx(1,1)/(gx+1.d-30))),    		     &
                      mesh%PGrid%dy(1,1)/(TraPar%uvp(TraPar%np)%v+1.d-30+    &
                      dsqrt(mesh%PGrid%dy(1,1)/(gy+1.d-30)))) 
     call TraPar%TrackingParticles(mesh%PGrid,mesh%PCell,Var,dt)
     time=time+dt
   end do
   romix=(row*mesh%PCell%vof(1,1)+roa*(1.d0-mesh%PCell%vof(1,1)))
   if(romix>1.5d0*TraPar%rop) then
     if(TraPar%Posp(TraPar%np)%x>pointp%x) call exit(1)
   else
     if(TraPar%Posp(TraPar%np)%x<pointp%x) call exit(1)
   end if  
END PROGRAM particletest

