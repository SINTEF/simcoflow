!FIXME: This has too many dependencies for now
Program clsvoftest
   USE PrecisionVar
   USE StateVariables, ONLY :TVariables,TWave
   USE Mesh, ONLY : TsimcoMesh, TPoint
   USE Clsvof, ONLY : SolidObject,frac,Coupled_LS_VOF
   USE Constants, ONLY : g,pi
   USE CutCell
   IMPLICIT NONE
   TYPE(TVariables) :: Var
   TYPE(TPoint)     :: point, point2, pointP, ReS, ReE
   TYPE(TsimcoMesh) :: mesh
   TYPE(TWave)      :: wave
   TYPE(SolidObject):: BoomCase
   INTEGER(it4b)    :: i, j, ii, jj
   REAL(dp)         :: dt, gx, gy, time, romix, gmix, Lw
   ! Need to create a mesh first to set Isize and Jsize used all over
   mesh = TsimcoMesh(10,12)
   Var = Tvariables()
   Var%p(:,:)=0.d0
   Var%u(:,:)=0.d0
   Var%v(:,:)=1.d0
   
   point = TPoint(0.d0,0.d0)
   point2 = TPoint(10.d0,10.d0)
   
   ReS = TPoint(2.d0,2.d0)
   ReE = TPoint(6.d0,6.d0)
   gx=0.d0
   gy=g
   
   
   call mesh%Initial_Grid2(point,point2,ReS,ReE,11,13,1,2,1.d0,0)
   call mesh%InitialUVGrid2(1.d0,1.d0)
   ! Initialize clsvof for PCell, UCell, VCell 
   Lw = 5.d0
   do i=1,sizeof(mesh%PGrid%x(:,1))
     do j=1,sizeof(mesh%PGrid%y(1,:))
       mesh%PCell%vofS(i,j)=0.d0
       mesh%PCell%phiS(i,j)=1.d4
       mesh%PCell%nxS(i,j)=1.d0
       mesh%PCell%nyS(i,j)=0.d0
       
       mesh%PCell%phi(i,j)=mesh%PGrid%y(i,j)-Lw
       if(mesh%PCell%phi(i,j)<Lw) then
         mesh%PCell%ny(i,j)=1.d0
       else
         mesh%PCell%ny(i,j)=-1.d0
       end if
       mesh%PCell%nx(i,j)=0.d0
       call frac(mesh%PCell%nx(i,j),mesh%PCell%ny(i,j),mesh%PCell%phi(i,j),    &
                 mesh%PGrid%dx(i,j),mesh%PGrid%dy(i,j),mesh%PCell%vof(i,j))

       mesh%UCell%vofS(i,j)=0.d0
       mesh%UCell%phiS(i,j)=1.d4
       mesh%UCell%nxS(i,j)=1.d0
       mesh%UCell%nyS(i,j)=0.d0
       
       mesh%UCell%phi(i,j)=mesh%UGrid%y(i,j)-Lw
       if(mesh%UCell%phi(i,j)<Lw) then
         mesh%UCell%ny(i,j)=1.d0
       else
         mesh%UCell%ny(i,j)=-1.d0
       end if
       mesh%UCell%nx(i,j)=0.d0
       call frac(mesh%UCell%nx(i,j),mesh%UCell%ny(i,j),mesh%UCell%phi(i,j),    &
                 mesh%UGrid%dx(i,j),mesh%UGrid%dy(i,j),mesh%UCell%vof(i,j))

       mesh%VCell%vofS(i,j)=0.d0
       mesh%VCell%phiS(i,j)=1.d4
       mesh%VCell%nxS(i,j)=1.d0
       mesh%VCell%nyS(i,j)=0.d0
       
       mesh%VCell%phi(i,j)=mesh%VGrid%y(i,j)-Lw
       if(mesh%VCell%phi(i,j)<Lw) then
         mesh%VCell%ny(i,j)=1.d0
       else
         mesh%VCell%ny(i,j)=-1.d0
       end if
       mesh%VCell%nx(i,j)=0.d0
       call frac(mesh%VCell%nx(i,j),mesh%VCell%ny(i,j),mesh%VCell%phi(i,j),    &
                 mesh%VGrid%dx(i,j),mesh%VGrid%dy(i,j),mesh%VCell%vof(i,j))
     end do
   end do
   ! Initialize the boom 
   BoomCase%Posp%x=5000.d0
   BoomCase%Posp%y=1000.d0
   BoomCase%Dobj=0.8d0
   BoomCase%Wobj=0.16d0
   BoomCase%XBar1=BoomCase%Posp%x-BoomCase%Wobj/2.d0
   BoomCase%XBar2=BoomCase%Posp%x+BoomCase%Wobj/2.d0
   BoomCase%LBar=1.5d0
   BoomCase%YBar=BoomCase%Posp%y-dsqrt((BoomCase%Dobj/2.d0)**2.d0-             &
                 (BoomCase%Wobj/2.d0)**2.d0)-BoomCase%LBar
   BoomCase%Mobj=(pi/4.d0*(BoomCase%Dobj)**2.d0+BoomCase%Wobj*BoomCase%LBar)*0.5d0  
   BoomCase%us=0.d0
   BoomCase%vs=0.d0
   ! 
   wave = TWave(1.d0,1.d0,0.d0,10.d0,2.d0,1.d0,20.d0,20.d0,2.d0,1.d0,5.d0)
   call Grid_Preprocess(mesh,Var,int8(1))    
   time = 0.d0
   do i=1,10
     dt=0.3d0*mesh%PGrid%dy(1,1)/Var%v(1,1) 
     call Coupled_LS_VOF(mesh%PGrid,mesh%PCell,mesh%UCell,mesh%VCell,Var,wave, &
                                               BoomCase,time,dt,INT8(i))
     time=time+dt
   end do
END PROGRAM clsvoftest

