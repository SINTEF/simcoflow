Module Clsvof
    USE PrecisionVar
    USE Mesh
    USE Cutcell
    USE StateVariables
    USE Matrix
    IMPLICIT NONE
    PRIVATE
    INTEGER,PARAMETER:: band_width = 4,nv = 10,nl = 5
    REAL(KIND=dp),PUBLIC,PARAMETER:: vofeps=1.d-14,tolp=1.d-10,TolDeno=1.d-24
    TYPE,PUBLIC:: SolidObject
      TYPE(Point):: Posp,PospO
      REAL(KIND=dp):: us,vs,asx,asy
      REAL(KIND=dp):: Dobj,Wobj,Mobj
      REAL(KIND=dp):: Xbar1,Xbar2,Ybar,Lbar
    END TYPE
    REAL(KIND=dp),DIMENSION(:,:),pointer:: vfl,vflS,nxs,nys
    REAL(KIND=dp),DIMENSION(:,:),pointer:: phi,phiS
    REAL(KIND=dp),DIMENSION(:,:),allocatable:: nx,ny
    REAL(KIND=dp),DIMENSION(8,2):: Vset
    PUBLIC:: Initial_Clsvof,Initial_ClsvofUV,Coupled_LS_VOF,ComputeForceObject,&
             SolidVolumeFraction
    Interface Initial_Clsvof
      Module procedure Initial_Clsvof
    End interface
    Interface Initial_ClsvofUV
      Module procedure Initial_ClsvofUV
    End interface
    Interface Coupled_LS_VOF
      Module procedure Coupled_LS_VOF
    End interface
    Interface ComputeForceObject
      Module procedure ComputeForceObject
    End interface
    Interface SolidVolumeFraction
      Module procedure SolidVolumeFraction
    End interface
    Contains
    SUBROUTINE Initial_Clsvof(TGrid,TCell,BoomCase)
      TYPE(Grid),INTENT(IN):: TGrid
      TYPE(Cell),INTENT(INOUT),target:: TCell
      TYPE(SolidObject),INTENT(IN):: BoomCase
      INTEGER(kind=it4b):: i,j
      REAL(KIND=dp):: dx,dy,dis,vol
      REAL(KIND=dp):: dx2,dy2,tol,fx,dfx
      REAL(KIND=dp),DIMENSION(:,:):: node(6,2),CutP(2,2),dpt(4)
      INTEGER(kind=it4b):: temp,templ,k
      REAL(KIND=dp):: epsil,nxx,nyy,nxy,vos,CylBar
      allocate(nx(Isize,Jsize))
      allocate(ny(Isize,Jsize))
   !  for wave only
      vfl => TCell%vof
      phi => TCell%phi
      phiS => TCell%phiS
      vflS => TCell%VofS
      epsil=1.d-24
      do j=1,Jsize
        do i=1,Isize
       !  for stay still wave as initial condition
       !   dis=TGrid%y(i,j)-Depthw
       !   nxx=0.d0
       !   nyy=1.d0
       !   call frac(nxx,nyy,dis,TGrid%dx(i,j),TGrid%dy(i,j),vol)
       ! *********
       !   vfl(i,j)=vol
       !   phi(i,j)=dis
       !   TCell%nx(i,j)=nxx
       !   TCell%ny(i,j)=nyy
       !  for sinusoidal wave
          temp=1
          templ=1
          dx2=TGrid%dx(i,j)/2.d0
          dy2=TGrid%dy(i,j)/2.d0
          dpt(1)=TGrid%y(i,j)-dy2-Depthw-Amp0*dsin(kw*(TGrid%x(i,j)+dx2))
          dpt(2)=TGrid%y(i,j)-dy2-Depthw-Amp0*dsin(kw*(TGrid%x(i,j)-dx2))
          dpt(3)=TGrid%y(i,j)+dy2-Depthw-Amp0*dsin(kw*(TGrid%x(i,j)-dx2))
          dpt(4)=TGrid%y(i,j)+dy2-Depthw-Amp0*dsin(kw*(TGrid%x(i,j)+dx2))
          if(dpt(1)>=0.d0) then
          node(temp,1)=dx2
          node(temp,2)=-dy2
          temp=temp+1
          end if
          if(dpt(1)*dpt(2)<0.d0) then
            node(temp,1)=TGrid%x(i,j)-dx2+TGrid%dx(i,j)*dabs(dpt(2))/          &
                                          (dabs(dpt(2))+dabs(dpt(1)))
            tol=1.d0
            do while(tol>1.d-13)
              fx=TGrid%y(i,j)-dy2-(Depthw+Amp0*dsin(kw*(node(temp,1))))
              dfx=-Amp0*kw*dcos(kw*(node(temp,1)))
              tol=dabs(fx/dfx)
              node(temp,1)=node(temp,1)-fx/dfx
            end do
            node(temp,1)=node(temp,1)-TGrid%x(i,j)
            node(temp,2)=-dy2
            CutP(templ,1)=node(temp,1)
            CutP(templ,2)=node(temp,2)
            templ=templ+1
            temp=temp+1
          end if
          if(dpt(2)>=0.d0) then
            node(temp,1)=-dx2
            node(temp,2)=-dy2
            temp=temp+1
          end if
          if(dpt(2)*dpt(3)<0.d0) then
            node(temp,1)=-dx2
            node(temp,2)=Depthw+Amp0*dsin(kw*(node(temp,1)+TGrid%x(i,j)))-     &
                                                           TGrid%y(i,j)
            CutP(templ,1)=node(temp,1)
            CutP(templ,2)=node(temp,2)
            templ=templ+1
            temp=temp+1
          end if
          if(dpt(3)>=0.d0) then
            node(temp,1)=-dx2
            node(temp,2)=dy2
            temp=temp+1
          end if
          if(dpt(3)*dpt(4)<0.d0) then
            node(temp,1)=TGrid%x(i,j)-dx2+TGrid%dx(i,j)*dabs(dpt(3))/          &
                                         (dabs(dpt(3))+dabs(dpt(4)))
            tol=1.d0
            do while(tol>1.d-13)
              fx=TGrid%y(i,j)+dy2-(Depthw+Amp0*dsin(kw*(node(temp,1))))
              dfx=-Amp0*kw*dcos(kw*(node(temp,1)))
              tol=dabs(fx/dfx)
              node(temp,1)=node(temp,1)-fx/dfx
            end do
            node(temp,1)=node(temp,1)-TGrid%x(i,j)
            node(temp,2)=dy2
            CutP(templ,1)=node(temp,1)
            CutP(templ,2)=node(temp,2)
            templ=templ+1
            temp=temp+1
          end if
          if(dpt(4)>=0.d0) then
            node(temp,1)=dx2
            node(temp,2)=dy2
            temp=temp+1
          end if
          if(dpt(4)*dpt(1)<0.d0) then
            node(temp,1)=dx2
            node(temp,2)=Depthw+Amp0*dsin(kw*(node(temp,1)+TGrid%x(i,j)))-     &
                                                           TGrid%y(i,j)
            CutP(templ,1)=node(temp,1)
            CutP(templ,2)=node(temp,2)
            templ=templ+1
            temp=temp+1
          end if
          node(temp,1)=node(1,1)
          node(temp,2)=node(1,2)
          vol=0.d0
          do k=1,temp-1
            vol=vol+0.5d0*(node(k,1)*node(k+1,2)-node(k+1,1)*node(k,2))
          end do
          vol=1.d0-dabs(vol/(TGrid%dx(i,j)*TGrid%dy(i,j)))
          if(templ==3) then
            nxx=CutP(2,2)-CutP(1,2)
            nyy=CutP(2,1)-CutP(1,1)
            vos=-Amp0*cw0*kw*dcos(kw*(TGrid%x(i,j)+TGrid%dx(i,j)))
            nxy=dsqrt(nxx**2.d0+nyy**2.d0)
            if(vos>0.d0) then
              nxx=-nxx/nxy
            else
              nxx=nxx/nxy
            end if
            nyy=dabs(nyy/nxy)
            dis=dabs(CutP(2,1)*CutP(1,2)-CutP(1,1)*CutP(2,2))/nxy*             &
                                                   dsign(1.d0,0.5d0-vol)
          else
            dis=TGrid%y(i,j)-(Depthw+Amp0*dsin(kw*(TGrid%x(i,j))))
            nxx=0.d0
            nyy=1.d0
          end if
          vfl(i,j)=vol
          phi(i,j)=dis
          TCell%nx(i,j)=nxx
          TCell%ny(i,j)=nyy
        end do
      end do
      CylBar=dsqrt((BoomCase%Dobj/2.d0)**2.d0-(BoomCase%Wobj/2.d0)**2.d0)
      do i=1,Isize
        do j=1,Jsize
     !  For boom cylinder
          dx=TGrid%x(i,j)-BoomCase%Posp%x
          dy=TGrid%y(i,j)-BoomCase%Posp%y
          phiS(i,j)=(dsqrt(dx**2.d0+dy**2.d0)-BoomCase%Dobj/2.d0)
          nx(i,j) = dx/dsqrt(dx**2.d0+dy**2.d0)
          ny(i,j) = dy/dsqrt(dx**2.d0+dy**2.d0)
     !  For region at left side of boom
          if(TGrid%x(i,j)<=BoomCase%XBar1) then
     !  For region under boom
            if(TGrid%y(i,j)<BoomCase%YBar) then
              dx=TGrid%x(i,j)-BoomCase%XBar1
              dy=TGRid%y(i,j)-BoomCase%YBar
              PhiS(i,j)=dsqrt(dx**2.d0+dy**2.d0)
              nx(i,j) = dx/dsqrt(dx**2.d0+dy**2.d0)
              ny(i,j) = dy/dsqrt(dx**2.d0+dy**2.d0)
     !  For region have same altitude as Boom bar
            else if(TGrid%y(i,j)<BoomCase%YBar+BoomCase%LBar) then
              dis=BoomCase%XBar1-TGrid%x(i,j)
              if(dabs(dis)<dabs(phiS(i,j))) then
                nx(i,j)=-1.d0
                ny(i,j)=0.d0
                phiS(i,j)=dis
              end if
            end if
      !  For region at right side of boom
          else if(TGrid%x(i,j)>=BoomCase%XBar2) then
      !  For region under boom
            if(TGrid%y(i,j)<BoomCase%YBar) then
              dx=TGrid%x(i,j)-BoomCase%XBar2
              dy=TGRid%y(i,j)-BoomCase%YBar
              PhiS(i,j)=dsqrt(dx**2.d0+dy**2.d0)
              nx(i,j)=dx/dsqrt(dx**2.d0+dy**2.d0)
              ny(i,j)=dy/dsqrt(dx**2.d0+dy**2.d0)
      !  For region have same altitude as Boom bar
            else if(TGrid%y(i,j)<BoomCase%YBar+BoomCase%LBar) then
              dis=TGrid%x(i,j)-BoomCase%XBar2
              if(dabs(dis)<dabs(phiS(i,j))) then
                nx(i,j)=1.d0
                ny(i,j)=0.d0
                phiS(i,j)=dis
              end if
            end if
       !  For region in side boom
          else
       !  For region under boom
            if(TGrid%y(i,j)<=BoomCase%YBar) then
              PhiS(i,j)=BoomCase%YBar-TGrid%y(i,j)
              nx(i,j)=0.d0
              ny(i,j)=-1.d0
       !  For region inside boom bar
            elseif(TGrid%y(i,j)<BoomCase%YBar+BoomCase%LBar) then
              dpt(1)=BoomCase%YBar-TGrid%y(i,j)
              dpt(2)=BoomCase%XBar1-TGrid%x(i,j)
              dpt(3)=TGrid%x(i,j)-BoomCase%XBar2
              if(dabs(dpt(1))<=dabs(dpt(2)).and.dabs(dpt(1))<=dabs(dpt(3))) then
                PhiS(i,j)=dpt(1)
                nx(i,j)=0.d0
                ny(i,j)=-1.d0
              end if
              if(dabs(dpt(2))<=dabs(dpt(1)).and.dabs(dpt(2))<=dabs(dpt(3))) then
                PhiS(i,j)=dpt(2)
                nx(i,j)=-1.d0
                ny(i,j)=0.d0
              end if
              if(dabs(dpt(3))<=dabs(dpt(1)).and.dabs(dpt(3))<=dabs(dpt(2))) then
                PhiS(i,j)=dpt(3)
                nx(i,j)=1.d0
                ny(i,j)=0.d0
              end if
         !  For region inside cylinder
            elseif(TGrid%y(i,j)<BoomCase%Posp%y) then
              if(dabs(TGrid%x(i,j)-0.5d0*(BoomCase%XBar1+BoomCase%XBar2))/     &
                dabs(TGrid%y(i,j)-BoomCase%Posp%y)<BoomCase%Wobj/2.d0/CylBar) then
                dpt(1)=-dsqrt((TGrid%x(i,j)-BoomCase%XBar1)**2.d0+             &
                              (TGrid%y(i,j)-BoomCase%YBar-BoomCase%LBar)**2.d0)
                dpt(2)=-dsqrt((TGrid%x(i,j)-BoomCase%XBar2)**2.d0+             &
                              (TGrid%y(i,j)-BoomCase%YBar-BoomCase%LBar)**2.d0)
                if(dabs(dpt(1))<dabs(dpt(2))) then
                  PhiS(i,j)=dpt(1)
                  nx(i,j)=(BoomCase%XBar1-TGrid%x(i,j))/dabs(dpt(1))
                  ny(i,j)=-(TGrid%y(i,j)-BoomCase%YBar-BoomCase%LBar)/dabs(dpt(1))
                else
                  PhiS(i,j)=dpt(2)
                  nx(i,j)=(BoomCase%XBar2-TGrid%x(i,j))/dabs(dpt(2))
                  ny(i,j)=-(TGrid%y(i,j)-BoomCase%YBar-BoomCase%LBar)/dabs(dpt(2))
                end if
              end if
            end if
          end if
          call frac(nx(i,j),ny(i,j),phiS(i,j),TGrid%dx(i,j),TGrid%dy(i,j),vol)
          vflS(i,j)=vol
       !  phiS(i,j)=1.d3
          if(vflS(i,j)<vofeps) vflS(i,j)=0.d0
          if(vflS(i,j)>=1.d0-vofeps) vflS(i,j)=1.d0

       !  The specific case for the bottom of boom
          if(vflS(i,j)>vofeps.and.vflS(i,j)<1.d0-vofeps) then
            if(TGrid%y(i,j)+TGrid%dy(i,j)/2.d0-tolp>BoomCase%YBar.and.         &
               TGrid%y(i,j)-TGrid%dy(i,j)/2.d0+tolp<BoomCase%YBar) then
              if(TGrid%x(i,j)-TGrid%dx(i,j)/2.d0+tolp<BoomCase%XBar1.and.      &
                 TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-tolp>BoomCase%XBar1) then
                CutP(1,1)=BoomCase%XBar1-TGrid%x(i,j)
                CutP(1,2)=TGrid%dy(i,j)/2.d0
                CutP(2,1)=TGrid%dx(i,j)/2.d0
                CutP(2,2)=BoomCase%YBar-TGrid%y(i,j)
                nxx=CutP(2,2)-CutP(1,2)
                nyy=CutP(2,1)-CutP(1,1)
                vflS(i,j)=0.5d0*dabs(nxx*nyy)/TGrid%dx(i,j)/TGrid%dy(i,j)
                nxy=dsqrt(nxx**2.d0+nyy**2.d0)
                nx(i,j)=-dabs(nxx/nxy)
                ny(i,j)=-dabs(nyy/nxy)
                phiS(i,j)=dabs(CutP(2,1)*CutP(1,2)-CutP(1,1)*CutP(2,2))/nxy
              end if
              if(TGrid%x(i,j)-TGrid%dx(i,j)/2.d0+tolp<BoomCase%XBar2.and.      &
                 TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-tolp>BoomCase%XBar2)then
                CutP(1,1)=BoomCase%XBar2-TGrid%x(i,j)
                CutP(1,2)=TGrid%dy(i,j)/2.d0
                CutP(2,1)=-TGrid%dx(i,j)/2.d0
                CutP(2,2)=BoomCase%YBar-TGrid%y(i,j)
                nxx=CutP(2,2)-CutP(1,2)
                nyy=CutP(2,1)-CutP(1,1)
                vflS(i,j)=0.5d0*dabs(nxx*nyy)/TGrid%dx(i,j)/TGrid%dy(i,j)
                nxy=dsqrt(nxx**2.d0+nyy**2.d0)
                nx(i,j)=dabs(nxx/nxy)
                ny(i,j)=-dabs(nyy/nxy)
                phiS(i,j)=dabs(CutP(2,1)*CutP(1,2)-CutP(1,1)*CutP(2,2))/nxy
              end if
            end if
          end if
       !  The specific area between boom cylinder and boom bar
          if(vflS(i,j)>vofeps.and.vflS(i,j)<1.d0-vofeps) then
            if(TGrid%y(i,j)+TGrid%dy(i,j)/2.d0-tolp>BoomCase%Posp%y-CylBar.and. &
               TGrid%y(i,j)-TGrid%dy(i,j)/2.d0+tolp<BoomCase%Posp%y-CylBar)then
            ! For the left side of object
              if(TGrid%x(i,j)-TGrid%dx(i,j)/2.d0+tolp<BoomCase%XBar1.and.       &
                 TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-tolp>BoomCase%XBar1) then
                CutP(1,1)=BoomCase%XBar1-TGrid%x(i,j)
                CutP(1,2)=-TGrid%dy(i,j)/2.d0
                if(dsqrt((TGrid%x(i,j)-TGrid%dx(i,j)/2.d0-BoomCase%Posp%x)**2.d0+&
                  (TGrid%y(i,j)+TGrid%dy(i,j)/2.d0-BoomCase%Posp%y)**2.d0)<     &
                  BoomCase%Dobj/2.d0) then
                  CutP(2,1)=-TGrid%dx(i,j)/2.d0
                  CutP(2,2)=dsqrt((BoomCase%Dobj/2.d0)**2.d0-                  &
                    (BoomCase%Wobj/2.d0+TGrid%dx(i,j)/2.d0)**2.d0)-            &
                    (BoomCase%Posp%y-TGrid%y(i,j))
                  nxx=CutP(2,2)-CutP(1,2)
                  nyy=CutP(2,1)-CutP(1,1)
                  vflS(i,j)=1.d0-0.5d0*dabs(nxx*nyy)/TGrid%dx(i,j)/TGrid%dy(i,j)
                  nxy=dsqrt(nxx**2.d0+nyy**2.d0)
                  nx(i,j)=-dabs(nxx/nxy)
                  ny(i,j)=-dabs(nyy/nxy)
                  phiS(i,j)=dabs(CutP(2,1)*CutP(1,2)-CutP(1,1)*CutP(2,2))/nxy* &
                            dsign(1.d0,0.5d0-vflS(i,j))
                else
                  CutP(2,1)=dsqrt((BoomCase%Dobj/2.d0)**2.d0-                  &
                    (BoomCase%Posp%y-TGrid%y(i,j)-TGrid%dy(i,j)/2.d0)**2.d0)-  &
                    (BoomCase%Posp%x-TGrid%x(i,j))
                  CutP(2,2)=TGrid%dy(i,j)/2.d0
                  nxx=CutP(2,2)-CutP(1,2)
                  nyy=CutP(2,1)-CutP(1,1)
                  vflS(i,j)=0.5d0*dabs(TGrid%dx(i,j)/2.d0-CutP(1,1)+           &
                                  TGrid%dx(i,j)/2.d0-CutP(2,1))/TGrid%dx(i,j)
                  nxy=dsqrt(nxx**2.d0+nyy**2.d0)
                  nx(i,j)=-dabs(nxx/nxy)
                  ny(i,j)=-dabs(nyy/nxy)
                  phiS(i,j)=dabs(CutP(2,1)*CutP(1,2)-CutP(1,1)*CutP(2,2))/nxy* &
                            dsign(1.d0,0.5d0-vflS(i,j))
                end if
              end if
            ! For the right side of object
              if(TGrid%x(i,j)-TGrid%dx(i,j)/2.d0+tolp<BoomCase%XBar2.and.      &
                 TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-tolp>BoomCase%XBar2)then
                CutP(1,1)=BoomCase%XBar2-TGrid%x(i,j)
                CutP(1,2)=-TGrid%dy(i,j)/2.d0
                if(dsqrt((TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-BoomCase%Posp%x)**2.d0+&
                  (TGrid%y(i,j)+TGrid%dy(i,j)/2.d0-BoomCase%Posp%y)**2.d0)<     &
                  BoomCase%Dobj/2.d0) then
                  CutP(2,1)=TGrid%dx(i,j)/2.d0
                  CutP(2,2)=dsqrt((BoomCase%Dobj/2.d0)**2.d0-                  &
                    (BoomCase%Wobj/2.d0+TGrid%dx(i,j)/2.d0)**2.d0)-            &
                    (BoomCase%Posp%y-TGrid%y(i,j))
                  nxx=CutP(2,2)-CutP(1,2)
                  nyy=CutP(2,1)-CutP(1,1)
                  vflS(i,j)=1.d0-0.5d0*dabs(nxx*nyy)/TGrid%dx(i,j)/TGrid%dy(i,j)
                  nxy=dsqrt(nxx**2.d0+nyy**2.d0)
                  nx(i,j)=dabs(nxx/nxy)
                  ny(i,j)=-dabs(nyy/nxy)
                  phiS(i,j)=dabs(CutP(2,1)*CutP(1,2)-CutP(1,1)*CutP(2,2))/nxy* &
                            dsign(1.d0,0.5d0-vflS(i,j))
                else
                  CutP(2,1)=dsqrt((BoomCase%Dobj/2.d0)**2.d0-                  &
                    (BoomCase%Posp%y-TGrid%y(i,j)-TGrid%dy(i,j)/2.d0)**2.d0)-  &
                    (BoomCase%Posp%x-TGrid%x(i,j))
                  CutP(2,2)=TGrid%dy(i,j)/2.d0
                  nxx=CutP(2,2)-CutP(1,2)
                  nyy=CutP(2,1)-CutP(1,1)
                  vflS(i,j)=0.5d0*dabs(CutP(1,1)+TGrid%dx(i,j)+CutP(2,1))/     &
                                       TGrid%dx(i,j)
                  nxy=dsqrt(nxx**2.d0+nyy**2.d0)
                  nx(i,j)=dabs(nxx/nxy)
                  ny(i,j)=-dabs(nyy/nxy)
                  phiS(i,j)=dabs(CutP(2,1)*CutP(1,2)-CutP(1,1)*CutP(2,2))/nxy* &
                            dsign(1.d0,0.5d0-vflS(i,j))
                end if
              end if
            end if
          end if
          TCell%nxS(i,j)=nx(i,j)
          TCell%nyS(i,j)=ny(i,j)
          vfl(i,j)=vfl(i,j)*(1.d0-vflS(i,j))
        end do
      end do
      nullify(vfl)
      nullify(phi)
      deallocate(nx)
      deallocate(ny)
    end subroutine Initial_Clsvof

    subroutine Initial_ClsVofUV(PCell,PGrid,TCell,TGrid,VPp,SPp,VPuv,SPuv,     &
                                                                    BoomCase,uv)
       IMPLICIT NONE
       TYPE(Cell),INTENT(IN):: PCell
       TYPE(Cell),INTENT(INOUT):: TCell
       TYPE(Grid),INTENT(IN):: PGrid,TGrid
       TYPE(SolidObject),INTENT(IN):: BoomCase
       REAL(KIND=dp),DIMENSION(:,:,:),allocatable,INTENT(IN)::SPp
       REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(INOUT)::VPp,VPuv,SPuv
       INTEGER(kind=it4b),INTENT(IN):: uv
       INTEGER(kind=it4b):: i,j
       REAL(KIND=dp):: posi(4),valu(4),PosTar,temn,tol
       REAL(KIND=dp):: volfl,volfr,volsl,volsr,phifl,phifr,phisl,phisr
       tol = 1.d-14
       call SolidVolumeFraction(TGrid,TCell,BoomCase)
       if(uv==0) then
         do j = jbeg,jbeg+Jsize-1
           do i = ibeg,ibeg+Isize-2
             SPuv(i,j)=0.d0
             VPuv(i,j)=0.d0
             if(i>=3.and.i<=Isize-2) then
               posi(1)=TGrid%x(i,j)-PGrid%x(i-1,j)
               posi(2)=TGrid%x(i,j)-PGrid%x(i,j)
               posi(3)=TGrid%x(i,j)-PGrid%x(i+1,j)
               posi(4)=TGrid%x(i,j)-PGrid%x(i+2,j)
               PosTar = 0.d0
             ! for liquid part
               Valu(1)=PCell%phi(i-1,j);Valu(2)=PCell%phi(i,j)
               Valu(3)=PCell%phi(i+1,j);Valu(4)=PCell%phi(i+2,j)
               TCell%phi(i,j)=LagrangePoly(posi,valu,PosTar)
             ! for liquid part: normal vector
               Valu(1)=PCell%nx(i-1,j);Valu(2)=PCell%nx(i,j)
               Valu(3)=PCell%nx(i+1,j);Valu(4)=PCell%nx(i+2,j)
               TCell%nx(i,j)=LagrangePoly(posi,valu,PosTar)
             ! for liquid part: normal vector
               Valu(1)=PCell%ny(i-1,j);Valu(2)=PCell%ny(i,j)
               Valu(3)=PCell%ny(i+1,j);Valu(4)=PCell%ny(i+2,j)
               TCell%ny(i,j)=LagrangePoly(posi,valu,PosTar)
             else
             ! for liquid part
               TCell%phi(i,j) = 0.5d0*(PCell%phi(i+1,j)+PCell%phi(i,j))
               TCell%nx(i,j) = 0.5d0*(PCell%nx(i+1,j)+PCell%nx(i,j))
               TCell%ny(i,j) = 0.5d0*(PCell%ny(i+1,j)+PCell%ny(i,j))
             end if
             temn=dsqrt(TCell%nx(i,j)**2.d0+TCell%ny(i,j)**2.d0)+1.d-14
             TCell%nx(i,j)=TCell%nx(i,j)/temn
             TCell%ny(i,j)=TCell%ny(i,j)/temn

             phifl=PCell%phi(i,j)+PGrid%dx(i,j)/4.d0*PCell%nx(i,j)
          !   phisl=PCell%phiS(i,j)+PGrid%dx(i,j)/4.d0*PCell%nxS(i,j)
             phisl=TCell%phiS(i,j)
           ! for a cell with only fluid phase
             if(PCell%vof(i,j)>=1.d0-vofeps.and.PCell%vofS(i,j)<epsi) then
               volfl=1.d0
          !     volsl=0.d0
             end if
           ! for a cell with fluid and air(gas)
             if(PCell%vof(i,j)<1.d0-vofeps.and.PCell%vof(i,j)>=vofeps.and.     &
                                                    PCell%vofS(i,j)<epsi) then
               call frac(PCell%nx(i,j),PCell%ny(i,j),phifl,PGrid%dx(i,j)/2.d0, &
                                                          PGrid%dy(i,j),volfl)
           !    volsl=0.d0
             end if
           ! for a cell with only air(gas)
             if(PCell%vof(i,j)<vofeps.and.PCell%vofS(i,j)<epsi) then
               volfl=0.d0
           !    volsl=0.d0
             end if
           ! for a cell with only solid
             if(PCell%vof(i,j)<vofeps.and.PCell%vofS(i,j)>=1.d0-epsi) then
               volfl=0.d0
           !    volsl=1.d0
             end if
           ! for a cell with air and solid
             if(PCell%vof(i,j)<vofeps.and.PCell%vofS(i,j)<1.d0-vofeps.and.     &
                                                PCell%vofS(i,j)>=epsi) then
               call frac(PCell%nxS(i,j),PCell%nyS(i,j),phisl,                  &
                                       PGrid%dx(i,j)/2.d0,PGrid%dy(i,j),volsl)
               volfl=0.d0
             end if
             if(PCell%vof(i,j)>=vofeps.and.PCell%vof(i,j)<1.d0-vofeps.and.     &
                PCell%vofS(i,j)>=epsi.and.PCell%vofS(i,j)<1.d0-epsi) then
               call CellGeoCal(PCell%nxS(i,j),PCell%nyS(i,j),phisl,            &
                               PCell%nx(i,j),PCell%ny(i,j),phifl,              &
                               PGrid%dx(i,j)/2.d0,PGrid%dy(i,j),volsl,volfl)
             end if
             if(i<Isize) then
               phifr=PCell%phi(i+1,j)-PGrid%dx(i+1,j)/4.d0*PCell%nx(i+1,j)
           !    phisr=PCell%phiS(i+1,j)-PGrid%dx(i+1,j)/4.d0*PCell%nxS(i+1,j)
               phisr=TCell%phiS(i,j)
           ! for a cell with only fluid phase
               if(PCell%vof(i+1,j)>=1.d0-vofeps.and.PCell%vofS(i+1,j)<epsi) then
                 volfr=1.d0
           !      volsr=0.d0
               end if
           ! for a cell with fluid and air(gas)
               if(PCell%vof(i+1,j)<1.d0-vofeps.and.PCell%vof(i+1,j)>=vofeps.and. &
                                                 PCell%vofS(i+1,j)<epsi) then
                 call frac(PCell%nx(i+1,j),PCell%ny(i+1,j),phifr,              &
                                   PGrid%dx(i+1,j)/2.d0,PGrid%dy(i+1,j),volfr)
           !      volsr=0.d0
               end if
           ! for a cell with only air(gas)
               if(PCell%vof(i+1,j)<vofeps.and.PCell%vofS(i+1,j)<epsi) then
                 volfr=0.d0
           !      volsr=0.d0
               end if
           ! for a cell with only solid
               if(PCell%vof(i+1,j)<vofeps.and.PCell%vofS(i+1,j)>=1.d0-epsi) then
                 volfr=0.d0
           !      volsr=1.d0
               end if
           ! for a cell with air and solid
               if(PCell%vof(i+1,j)<vofeps.and.PCell%vofS(i+1,j)<1.d0-epsi.and. &
                                            PCell%vofS(i+1,j)>=epsi) then
                 call frac(PCell%nxS(i+1,j),PCell%nyS(i+1,j),phisr,            &
                                    PGrid%dx(i+1,j)/2.d0,PGrid%dy(i+1,j),volsr)
                 volfr=0.d0
               end if
           ! for a cell with air, fluid and solid
               if(PCell%vof(i+1,j)>=vofeps.and.PCell%vof(i+1,j)<1.d0-vofeps.and. &
                  PCell%vofS(i+1,j)>=epsi.and.PCell%vofS(i+1,j)<1.d0-epsi) then
                 call CellGeoCal(PCell%nxS(i+1,j),PCell%nyS(i+1,j),phisr,      &
                                PCell%nx(i+1,j),PCell%ny(i+1,j),phifr,         &
                                PGrid%dx(i+1,j)/2.d0,PGrid%dy(i+1,j),volsr,volfr)
               end if
             else
               volfr=volfl
          !     volsr=volsl
             end if
             TCell%vof(i,j)=0.5d0*(volfl+volfr)
             if(PCell%vof(i,j)>vofeps.and.SPp(i,j,1)>tol) then
               SPuv(i,j)=SPuv(i,j)+SPp(i,j,1)*volfl/PCell%vof(i,j)
               VPuv(i,j)=VPuv(i,j)+VPp(i,j)*volfl/PCell%vof(i,j)
             end if
             if(PCell%vof(i+1,j)>vofeps.and.SPp(i+1,j,1)>tol) then
               SPuv(i,j)=SPuv(i,j)+SPp(i+1,j,1)*volfr/PCell%vof(i+1,j)
               VPuv(i,j)=VPuv(i,j)+VPp(i+1,j)*volfr/PCell%vof(i+1,j)
             end if
             if(TCell%vof(i,j)<vofeps) TCell%vof(i,j)=0.d0
             if(TCell%vof(i,j)>1.d0-vofeps-TCell%vofS(i,j)) then
               TCell%vof(i,j)=1.d0-TCell%vofS(i,j)
             end if
             if(i>1) then
               if((PCell%vof(i-1,j)>=1.d0-vofeps-PCell%vofS(i-1,j).and.          &
                 PCell%vof(i+1,j)>=1.d0-vofeps-PCell%vofS(i+1,j))) then
                 TCell%vof(i,j)=1.d0-TCell%vofS(i,j)
               end if
             end if
         !    TCell%vofS(i,j)=0.5d0*(volsl+volsr)
         !    if(TCell%vofS(i,j)<epsi) TCell%VofS(i,j)=0.d0
         !    if(TCell%vofS(i,j)>=1.d0-epsi) TCell%VofS(i,j)=1.d0
           end do
           TCell%phi(i,j) = TCell%phi(i-1,j)
           TCell%nx(i,j) = TCell%nx(i-1,j)
           TCell%ny(i,j) = TCell%ny(i-1,j)
           TCell%vof(i,j) = TCell%vof(i-1,j)
         end do
       else
         do i = ibeg,ibeg+Isize-1
           do j = jbeg,jbeg+Jsize-2
             SPuv(i,j)=0.d0
             if(j>=3.and.j<=Jsize-2) then
               posi(1)=TGrid%y(i,j)-PGrid%y(i,j-1)
               posi(2)=TGrid%y(i,j)-PGrid%y(i,j)
               posi(3)=TGrid%y(i,j)-PGrid%y(i,j+1)
               posi(4)=TGrid%y(i,j)-PGrid%y(i,j+2)
               PosTar = 0.d0
             ! liquid level-set function for bottom half of cell
               Valu(1)=PCell%phi(i,j-1);Valu(2)=PCell%phi(i,j)
               Valu(3)=PCell%phi(i,j+1);Valu(4)=PCell%phi(i,j+2)
               TCell%phi(i,j)=LagrangePoly(posi,valu,PosTar)
             ! liquid normal vector for bottom half of cell
               Valu(1)=PCell%nx(i,j-1);Valu(2)=PCell%nx(i,j)
               Valu(3)=PCell%nx(i,j+1);Valu(4)=PCell%nx(i,j+2)
               TCell%nx(i,j)=LagrangePoly(posi,valu,PosTar)
             ! liquid normal vector for bottom half of cell
               Valu(1)=PCell%ny(i,j-1);Valu(2)=PCell%ny(i,j)
               Valu(3)=PCell%ny(i,j+1);Valu(4)=PCell%ny(i,j+2)
               TCell%ny(i,j)=LagrangePoly(posi,valu,PosTar)
             else
               TCell%phi(i,j)=0.5d0*(PCell%phi(i,j+1)+PCell%phi(i,j))
               TCell%nx(i,j)=0.5d0*(PCell%nx(i,j+1)+PCell%nx(i,j))
               TCell%ny(i,j)=0.5d0*(PCell%ny(i,j+1)+PCell%ny(i,j))
             end if
             temn=dsqrt(TCell%nx(i,j)**2.d0+TCell%ny(i,j)**2.d0+1.d-30)
             TCell%nx(i,j)=TCell%nx(i,j)/temn
             TCell%ny(i,j)=TCell%ny(i,j)/temn

             phifl=PCell%phi(i,j)+PGrid%dy(i,j)/4.d0*PCell%ny(i,j)
           !  phisl=PCell%phiS(i,j)+PGrid%dy(i,j)/4.d0*PCell%nyS(i,j)
             phisl=TCell%phiS(i,j)
           ! for a cell with only fluid phase
             if(PCell%vof(i,j)>=1.d0-vofeps.and.PCell%vofS(i,j)<epsi) then
               volfl=1.d0
           !    volsl=0.d0
             end if
           ! for a cell with fluid and air(gas)
             if(PCell%vof(i,j)<1.d0-vofeps.and.PCell%vof(i,j)>=vofeps.and.     &
                                                    PCell%vofS(i,j)<epsi) then
               call frac(PCell%nx(i,j),PCell%ny(i,j),phifl,PGrid%dx(i,j),      &
                                                    PGrid%dy(i,j)/2.d0,volfl)
           !    volsl=0.d0
             end if
           ! for a cell with only air(gas)
             if(PCell%vof(i,j)<vofeps.and.PCell%vofS(i,j)<epsi) then
               volfl=0.d0
           !    volsl=0.d0
             end if
           ! for a cell with only solid
             if(PCell%vof(i,j)<vofeps.and.PCell%vofS(i,j)>=1.d0-epsi) then
               volfl=0.d0
           !    volsl=1.d0
             end if
           ! for a cell with air and solid
             if(PCell%vof(i,j)<vofeps.and.PCell%vofS(i,j)<1.d0-epsi.and.       &
                                              PCell%vofS(i,j)>=epsi) then
               call frac(PCell%nxS(i,j),PCell%nyS(i,j),phisl,                  &
                                     PGrid%dx(i,j),PGrid%dy(i,j)/2.d0,volsl)
               volfl=0.d0
             end if
           ! for a cell with air, fluid, solid
             if(PCell%vof(i,j)>=vofeps.and.PCell%vof(i,j)<1.d0-vofeps.and.     &
                PCell%vofS(i,j)>=epsi.and.PCell%vofS(i,j)<1.d0-epsi) then
               call CellGeoCal(PCell%nxS(i,j),PCell%nyS(i,j),phisl,            &
                               PCell%nx(i,j),PCell%ny(i,j),phifl,              &
                               PGrid%dx(i,j),PGrid%dy(i,j)/2.d0,volsl,volfl)
             end if

             if(j<Jsize) then
               phifr=PCell%phi(i,j+1)-PGrid%dy(i,j+1)/4.d0*PCell%ny(i,j+1)
           !    phisr=PCell%phiS(i,j+1)-PGrid%dy(i,j+1)/4.d0*PCell%nyS(i,j+1)
               phisr=TCell%phiS(i,j)
           ! for a cell with only fluid phase
               if(PCell%vof(i,j+1)>=1.d0-vofeps.and.PCell%vofS(i,j+1)<epsi) then
                 volfr=1.d0
          !       volsr=0.d0
               end if
           ! for a cell with fluid and air(gas)
               if(PCell%vof(i,j+1)<1.d0-vofeps.and.PCell%vof(i,j+1)>=vofeps.and. &
                                                 PCell%vofS(i,j+1)<epsi) then
                 call frac(PCell%nx(i,j+1),PCell%ny(i,j+1),phifr,              &
                                   PGrid%dx(i,j+1),PGrid%dy(i,j+1)/2.d0,volfr)
           !      volsr=0.d0
               end if
           ! for a cell with only air(gas)
               if(PCell%vof(i,j+1)<vofeps.and.PCell%vofS(i,j+1)<epsi) then
                 volfr=0.d0
           !      volsr=0.d0
               end if
           ! for a cell with only solid
               if(PCell%vof(i,j+1)<vofeps.and.PCell%vofS(i,j+1)>=1.d0-epsi) then
                 volfr=0.d0
           !      volsr=1.d0
               end if
           ! for a cell with air and solid
               if(PCell%vof(i,j+1)<vofeps.and.PCell%vofS(i,j+1)<1.d0-epsi.and. &
                                            PCell%vofS(i,j+1)>=epsi) then
                 call frac(PCell%nxS(i,j+1),PCell%nyS(i,j+1),phisr,            &
                                    PGrid%dx(i,j+1),PGrid%dy(i,j+1)/2.d0,volsr)
                 volfr=0.d0
               end if
           ! for a cell with air, fluid and solid
               if(PCell%vof(i,j+1)>=vofeps.and.PCell%vof(i,j+1)<1.d0-vofeps.and. &
                  PCell%vofS(i,j+1)>=epsi.and.PCell%vofS(i,j+1)<1.d0-epsi) then
                 call CellGeoCal(PCell%nxS(i,j+1),PCell%nyS(i,j+1),phisr,      &
                                PCell%nx(i,j+1),PCell%ny(i,j+1),phifr,         &
                                PGrid%dx(i,j+1),PGrid%dy(i,j+1)/2.d0,volsr,volfr)
               end if
             else
               volfr=volfl
           !    volsr=volsl
             end if
             TCell%vof(i,j)=0.5d0*(volfl+volfr)
             if(PCell%vof(i,j)>vofeps.and.SPp(i,j,2)>tol) then
               SPuv(i,j)=SPuv(i,j)+SPp(i,j,2)*volfl/PCell%vof(i,j)
               VPuv(i,j)=VPuv(i,j)+VPp(i,j)*volfl/PCell%vof(i,j)
             end if
             if(PCell%vof(i,j+1)>vofeps.and.SPp(i,j+1,2)>tol) then
               SPuv(i,j)=SPuv(i,j)+SPp(i,j+1,2)*volfr/PCell%vof(i,j+1)
               VPuv(i,j)=VPuv(i,j)+VPp(i,j+1)*volfr/PCell%vof(i,j+1)
             end if
             if(TCell%vof(i,j)<vofeps) TCell%vof(i,j)=0.d0
             if(TCell%vof(i,j)>1.d0-vofeps-TCell%vofS(i,j))then
               TCell%vof(i,j)=1.d0-TCell%vofS(i,j)
             end if
             if(j>1) then
               if((PCell%vof(i,j-1)>=1.d0-vofeps-PCell%vofS(i,j-1).and.           &
                 PCell%vof(i,j+1)>=1.d0-vofeps-PCell%vofS(i,j+1))) then
                 TCell%vof(i,j)=1.d0-TCell%vofS(i,j)
               end if
             end if
          !   TCell%vofS(i,j)=0.5d0*(volsl+volsr)
          !   if(TCell%vofS(i,j)<epsi) TCell%vofS(i,j)=0.d0
         !    if(TCell%vofS(i,j)>=1.d0-epsi) TCell%VofS(i,j)=1.d0
           end do
         ! Boundary condition for liquid part
           TCell%phi(i,j)=TCell%phi(i,j-1)
           TCell%nx(i,j)=TCell%nx(i,j-1)
           TCell%ny(i,j)=TCell%ny(i,j-1)
           TCell%vof(i,j)=TCell%vof(i,j-1)
         ! Boundary condition for solid part
       !    TCell%phiS(i,j)=TCell%phiS(i,j-1)
       !    TCell%nxS(i,j)=TCell%nxS(i,j-1)
       !    TCell%nyS(i,j)=TCell%nyS(i,j-1)
       !    TCell%vofS(i,j)=TCell%vofS(i,j-1)
         end do
       end if
    end subroutine Initial_ClsVofUV

    function LagrangePoly(posi,valu,PosTar) result(ValTar)
       IMPLICIT NONE
       REAL(KIND=dp):: posi(4),valu(4),PosTar,ValTar,Mult
       INTEGER:: i,j
       ValTar = 0.d0
       do i = 1,4
         Mult = valu(i)
         do j = 1,4
           if(j/=i) then
             Mult = Mult*(PosTar-Posi(j))/(Posi(i)-Posi(j))
           end if
         end do
         ValTar = ValTar+Mult
       end do
    end function LagrangePoly
!***********************************************************
! Coupled level set and volume of fluid
!***********************************************************
    subroutine Coupled_LS_VOF(PGrid,PCell,UCell,VCell,TVar,BoomCase,NondiT,dt,itt)
       IMPLICIT NONE
       TYPE(Grid),INTENT(IN):: PGrid
       TYPE(Cell),INTENT(IN):: UCell,VCell
       TYPE(Cell),INTENT(INOUT),target:: PCell
       TYPE(Variables),INTENT(IN):: TVar
       TYPE(SolidObject),INTENT(INOUT):: BoomCase
       REAL(KIND=dp),INTENT(IN):: dt,NondiT
       INTEGER(kind=it8b),INTENT(IN):: itt
       INTEGER(kind=it4b):: i,j,ii,jj,k,nv,UpdateNorVec,tempc
       REAL(KIND=dp):: dtv,epsil
       REAL(KIND=dp),DIMENSION(:,:),allocatable:: dis,temvfx,temvfy,temlsx,temlsy
       REAL(KIND=dp),DIMENSION(:,:),allocatable:: ue,ve
       REAL(KIND=dp):: weight(9),sumvol,sumweight,sumweight2,neighvol,voldif
       REAL(KIND=dp):: SumAllVolBig,SumAllVolSmall,ResVol
       REAL(KIND=dp):: del_x,del_y,dttol,vb
       Vset(1,1)=0.d0;Vset(1,2)=-1.d0
       Vset(2,1)=0.d0;Vset(1,2)=1.d0
       Vset(3,1)=-1.d0;Vset(3,2)=0.d0
       Vset(4,1)=-1.d0;Vset(4,2)=-1.d0
       Vset(5,1)=-1.d0;Vset(5,2)=1.d0
       Vset(6,1)=1.d0;Vset(6,2)=0.d0
       Vset(7,1)=1.d0;Vset(7,2)=-1.d0
       Vset(8,1)=1.d0;Vset(8,2)=1.d0
       nv = 6
       vfl=>PCell%vof
       phi=>PCell%phi
       nxs=>PCell%nxs
       nys=>PCell%nys
       vflS=>PCell%VofS
       phiS=>PCell%phiS
       dtv=dt/dble(nv)
       dttol=dt !dtv
       del_x=PGrid%dx(Isize/2,Jsize/2)
       del_y=PGrid%dy(Isize/2,Jsize/2)
       epsil=dsqrt(2.d0)*del_x/2.d-6
       ResVol=1.d0
       allocate(nx(Isize,Jsize))
       allocate(ny(Isize,Jsize))
       allocate(dis(Isize,Jsize))
       allocate(ue(0:Isize,Jsize))
       allocate(ve(Isize,0:Jsize))
       allocate(temvfx(Isize,Jsize))
       allocate(temvfy(Isize,Jsize))
       allocate(temlsx(Isize,Jsize))
       allocate(temlsy(Isize,Jsize))
       BoomCase%PospO=BoomCase%Posp
       vb=BoomCase%vs
       do j = 1,JSize
         do i = 1,ISize
           if(UCell%MoExCell(i,j)/=1.and.UCell%Posnu(i,j)/=-1)then
             ue(i,j)=TVar%u(i,j)
           else
             ue(i,j)=0.d0
           end if
         end do
         ue(0,j)=TVar%u(0,j)
       end do
       do i = 1,ISize
         do j = 1,JSize
           if(VCell%MoExCell(i,j)/=1.and.VCell%Posnu(i,j)/=-1) then
             ve(i,j)=TVar%v(i,j)
           else
             ve(i,j)=0.d0
           end if
         end do
         ve(i,0)=TVar%v(i,0)
       end do
       do k = 1,nv
         if(mod(k,2)==1) then
           temvfx(:,:) = vfl(:,:)
           temlsx(:,:) = phi(:,:)
           UpdateNorVec=1
           call X_Sweep(PGrid,PCell,temvfx,temlsx,ue,ve,nx,ny,dis,dtv,itt,     &
                                                          UpdateNorVec,NondiT,k)
           do i = 1,ISize
             do j = 1,JSize
               temvfx(i,j)=temvfx(i,j)/(1.d0-dtv/PGrid%dx(i,j)*                &
                    (ue(i,j)*PCell%EEdge_Area(i,j)-ue(i-1,j)*PCell%WEdge_Area(i,j)))
               temlsx(i,j)=temlsx(i,j)/(1.d0-dtv/PGrid%dx(i,j)*                &
                    (ue(i,j)*PCell%EEdge_Area(i,j)-ue(i-1,j)*PCell%WEdge_Area(i,j)))
               vfl(i,j)=temvfx(i,j)*(1.d0-vflS(i,j))
               phi(i,j)=temlsx(i,j)
               if(vfl(i,j)>=1.d0-vofeps-vfls(i,j)) vfl(i,j)=1.d0-vfls(i,j)
               if(vfl(i,j)<vofeps) vfl(i,j) = 0.d0
               if(isnan(vfl(i,j))) then
                 print*,'what the fuck'
                 print*,temvfx(i,j),temvfy(i,j)
                 print*,'clsvof_823'
                 print*,i,j
               end if
           !    if(itt==3.and.i==1.and.j==149) then
           !      print*,itt,k
           !      print*,temvfx(i,j),temvfx(i,j)*(1.d0-dtv/PGrid%dx(i,j)*        &
           !         (ue(i,j)*PCell%EEdge_Area(i,j)-ue(i-1,j)*PCell%WEdge_Area(i,j)))
           !      print*,ue(i,j),ue(i-1,j)
           !      print*,(1.d0-dtv/PGrid%dx(i,j)*                                &
           !         (ue(i,j)*PCell%EEdge_Area(i,j)-ue(i-1,j)*PCell%WEdge_Area(i,j)))
           !      print*,vfl(i,j)
           !      print*,
          !     end if
             end do
           end do
       !    call Boundary_Condition_Vof_Phi(PGrid,NondiT+dble(k)*dtv)
           temvfy(:,:) = vfl(:,:)
           temlsy(:,:) = phi(:,:)
           UpdateNorVec=1
           call Y_Sweep(PGrid,PCell,temvfy,temlsy,ue,ve,nx,ny,dis,dtv,itt,     &
                                                                  UpdateNorVec)
           do i = 1,ISize
             do j = 1,JSize
               temvfy(i,j)=temvfy(i,j)+dtv/PGrid%dy(i,j)*temvfx(i,j)*          &
                (ve(i,j)*PCell%NEdge_Area(i,j)-ve(i,j-1)*PCell%SEdge_Area(i,j))!&
                ! -temvfx(i,j)*vb*dtv*PCell%nyS(i,j)*PCell%WlLh(i,j)/           &
                !  PGrid%dx(i,j)/PGrid%dy(i,j)
               temlsy(i,j)=temlsy(i,j)+dtv/PGrid%dy(i,j)*temlsx(i,j)*          &
                (ve(i,j)*PCell%NEdge_Area(i,j)-ve(i,j-1)*PCell%SEdge_Area(i,j))!&
                ! -temlsx(i,j)*vb*dtv*PCell%nyS(i,j)*PCell%WlLh(i,j)/           &
                !  PGrid%dx(i,j)/PGrid%dy(i,j)
          !    Correction for mass error conservation(Mark Sussmann JCP 221(2007) 469-505)
               vfl(i,j)=temvfy(i,j)*(1.d0-vflS(i,j))
               phi(i,j)=temlsy(i,j)
         !      vfl(i,j)=(temvfy(i,j)-temvfx(i,j)*dtv*                          &
         !!                          ((ue(i,j)*PCell%EEdge_Area(i,j)-            &
         !                            ue(i-1,j)*PCell%WEdge_Area(i,j))/         &
         !             PGrid%dx(i,j)+(ve(i,j)*PCell%NEdge_Area(i,j)-            &
         !             ve(i,j-1)*PCell%SEdge_Area(i,j))/PGrid%dy(i,j)))*(1.d0-vflS(i,j))
         !      phi(i,j)=temlsy(i,j)-temlsx(i,j)*dtv*                           &
         !                          ((ue(i,j)*PCell%EEdge_Area(i,j)-            &
         !                            ue(i-1,j)*PCell%WEdge_Area(i,j))/         &
         !             PGrid%dx(i,j)+(ve(i,j)*PCell%NEdge_Area(i,j)-            &
         !             ve(i,j-1)*PCell%SEdge_Area(i,j))/PGrid%dy(i,j))
          !     if(vfl(i,j)>=1.d0-vofeps-vfls(i,j)) vfl(i,j)=1.d0-vfls(i,j)
          !     if(vfl(i,j)<vofeps) vfl(i,j) = 0.d0
          !     if(itt==361.and.(i==2.or.i==3).and.(j==155.or.j==156.or.j==157.or.j==158.or.j==159)) then
         !      end if
            !   if(itt==3.and.i==1.and.j==149) then
            !     print*,itt,k
            !     print*,temvfy(i,j),temvfy(i,j)-dtv/PGrid%dy(i,j)*temvfx(i,j)*          &
            !    (ve(i,j)*PCell%NEdge_Area(i,j)-ve(i,j-1)*PCell%SEdge_Area(i,j))
            !     print*,vfl(i,j)
            !     print*,
            !   end if
               if(isnan(vfl(i,j))) then
                 print*,'what the fuck'
                 print*,temvfx(i,j),temvfy(i,j)
                 print*,'clsvof_858'
                 print*,i,j
               end if
             end do
           end do
         else
           temvfy(:,:)=vfl(:,:)
           temlsy(:,:)=phi(:,:)
           UpdateNorVec=1
           call Y_Sweep(PGrid,PCell,temvfy,temlsy,ue,ve,nx,ny,dis,dtv,itt,     &
                                                                   UpdateNorVec)
           do i = 1,ISize
             do j = 1,JSize
               temvfy(i,j)=temvfy(i,j)/(1.d0-dtv/PGrid%dy(i,j)*                &
                 (ve(i,j)*PCell%NEdge_Area(i,j)-ve(i,j-1)*PCell%SEdge_Area(i,j)))
               temlsy(i,j)=temlsy(i,j)/(1.d0-dtv/PGrid%dy(i,j)*                &
                 (ve(i,j)*PCell%NEdge_Area(i,j)-ve(i,j-1)*PCell%SEdge_Area(i,j)))
               vfl(i,j)=temvfy(i,j)*(1.d0-vflS(i,j))
               phi(i,j)=temlsy(i,j)
               if(vfl(i,j)>=1.d0-vofeps-vfls(i,j)) vfl(i,j)=1.d0-vfls(i,j)
               if(vfl(i,j)<vofeps) vfl(i,j)=0.d0
            !   if(itt==45.and.i==299.and.(j==123.or.j==124)) then
            !     print*,'test Y-Y'
            !     print*,k
            !     print*,i,j
            !     print*,temvfx(i,j),temvfy(i,j)
            !     print*,vfl(i,j)
            !     print*,ve(i,j),PCell%NEdge_Area(i,j)
            !     print*,ve(i,j-1),PCell%SEdge_Area(i,j)
            !    print*,
            !   end if
               if(isnan(vfl(i,j))) then
                 print*,'what the fuck'
                 print*,temvfx(i,j),temvfy(i,j)
                 print*,'clsvof_888'
                 print*,i,j
               end if
             end do
           end do
       !    call Boundary_Condition_Vof_Phi(PGrid,NondiT+dble(k)*dtv)
           temvfx(:,:) = vfl(:,:)
           temlsx(:,:) = phi(:,:)
           UpdateNorVec=1
           call X_Sweep(PGrid,PCell,temvfx,temlsx,ue,ve,nx,ny,dis,dtv,itt,     &
                                                         UpdateNorVec,NondiT,k)
           do i = 1,ISize
             do j = 1,JSize
               temvfx(i,j)=temvfx(i,j)+dtv/PGrid%dx(i,j)*temvfy(i,j)*          &
                (ue(i,j)*PCell%EEdge_Area(i,j)-ue(i-1,j)*PCell%WEdge_Area(i,j))!& ! By using the FVM the contribution of wall movement including in the transport equation
          !      -temvfy(i,j)*vb*dtv*PCell%nyS(i,j)*PCell%WlLh(i,j)/            &
          !       PGrid%dx(i,j)/PGrid%dy(i,j)
               temlsx(i,j)=temlsx(i,j)+dtv/PGrid%dx(i,j)*temlsy(i,j)*          &
                (ue(i,j)*PCell%EEdge_Area(i,j)-ue(i-1,j)*PCell%WEdge_Area(i,j))!&
         !       -temlsy(i,j)*vb*dtv*PCell%nyS(i,j)*PCell%WlLh(i,j)/            &
         !        PGrid%dx(i,j)/PGrid%dy(i,j)
         ! Additional step to maintain mass conservation
               vfl(i,j)=temvfx(i,j)*(1.d0-vflS(i,j))
               phi(i,j)=temlsx(i,j)
         !      vfl(i,j)=(temvfx(i,j)-temvfy(i,j)*dtv*                          &
         !                          ((ue(i,j)*PCell%EEdge_Area(i,j)-            &
         !                            ue(i-1,j)*PCell%WEdge_Area(i,j))/         &
         !             PGrid%dx(i,j)+(ve(i,j)*PCell%NEdge_Area(i,j)-            &
         !             ve(i,j-1)*PCell%SEdge_Area(i,j))/PGrid%dy(i,j)))*(1.d0-vflS(i,j))
         !      phi(i,j)=temlsx(i,j)-temlsx(i,j)*dtv*                           &
         !                          ((ue(i,j)*PCell%EEdge_Area(i,j)-            &
         !                            ue(i-1,j)*PCell%WEdge_Area(i,j))/         &
         !             PGrid%dx(i,j)+(ve(i,j)*PCell%NEdge_Area(i,j)-            &
         !             ve(i,j-1)*PCell%SEdge_Area(i,j))/PGrid%dy(i,j))
         !      if(vfl(i,j)>=1.d0-vofeps-vfls(i,j)) vfl(i,j)=1.d0-vfls(i,j)
         !      if(vfl(i,j)<vofeps) vfl(i,j) = 0.d0
         !      end if
               if(isnan(vfl(i,j))) then
                 print*,'what the fuck'
                 print*,temvfx(i,j),temvfy(i,j)
                 print*,'clsvof_929'
                 print*,i,j
               end if
             end do
           end do
         end if
         call ObjectMovement(BoomCase,dtv)
         call SolidVolumeFraction(Pgrid,PCell,BoomCase)
         call Boundary_Condition(vfl)
         call Boundary_Condition(phi)
       !  call Boundary_Condition_Vof_Phi(PGrid,NondiT+dble(k)*dtv)
       ! this step will update the level-set field by reconstructing interface
       ! the outputs are interface normal vector and a distance from cell center to interface
       ! Redistribute the volume of fluid to neighbor cells
       ! For j=1
         do i = 2,Isize-1
           if(vfl(i,1)>=1.d0-vfls(i,1)+vofeps) then
             sumvol=0.d0
             tempc=1
             weight=0.d0
             do ii=-1,1
               do jj=0,1
                 weight(tempc)=dmin1(ResVol,dmax1(0.d0,1.d0-vfls(i+ii,1+jj)-   &
                    vfl(i+ii,1+jj)))*(2.d0-dabs(dble(ii))-dabs(dble(jj))+      &
                    vfls(i+ii,1+jj)/(vfls(i+ii,1+jj)+1.d-15)*1.d-5)
                 tempc=tempc+1
               end do
             end do
             sumweight=sum(weight)+1.d-20-weight(3)
             weight(3)=0.d0
             voldif=(vfl(i,1)+vfls(i,1)-1.d0)
             tempc=1
             do ii=-1,1
               do jj=0,1
                 neighvol=voldif*weight(tempc)/sumweight
                 sumvol=sumvol+neighvol
                 vfl(i+ii,1+jj)=vfl(i+ii,1+jj)+neighvol
                 tempc=tempc+1
               end do
             end do
             vfl(i,1)=vfl(i,1)-sumvol
           end if
         end do

         do j = 2,Jsize-1
           do i = 2,Isize-1
             if(vfl(i,j)>=1.d0-vfls(i,j)+vofeps) then
               sumvol=0.d0
               tempc=1
               do ii=-1,1
                 do jj=-1,1
                   weight(tempc)=dmin1(ResVol,dmax1(0.d0,1.d0-vfls(i+ii,j+jj)- &
                     vfl(i+ii,j+jj)))*(2.d0-dabs(dble(ii))-dabs(dble(jj))+     &
                     vfls(i+ii,j+jj)/(vfls(i+ii,j+jj)+1.d-15)*1.d-5)
                   tempc=tempc+1
                 end do
               end do
             ! Flux out to neighbor cell
               sumweight=sum(weight)+1.d-20-weight(5)
               weight(5)=0.d0
               voldif=(vfl(i,j)+vfls(i,j)-1.d0)
               tempc=1
               do ii=-1,1
                 do jj=-1,1
                   neighvol=voldif*weight(tempc)/sumweight
                   sumvol=sumvol+neighvol
                   vfl(i+ii,j+jj)=vfl(i+ii,j+jj)+neighvol
                   tempc=tempc+1
                 end do
               end do
               vfl(i,j)=vfl(i,j)-sumvol
             end if
           end do
         end do

        ! Calculate the residual of vfl in cells that excess 1.d0-vfls
        ! Calculate the deficit of vfl in cells that is smaller than 1.d0-vfls.
        ! however, their value should be 1.d0-vfls
         sumvol=0.d0
         sumweight=0.d0
         tempc=1
         do j = 1,Jsize
           do i = 1,Isize
             if(vfl(i,j)>=1.d0-vfls(i,j)+vofeps) then
            !   if(i==165.and.j==160) then
            !     print*,vfl(i,j),1.d0-vfls(i,j)+vofeps
            !     print*,'before test'
            !   end if
               sumvol=sumvol+(vfl(i,j)+vfls(i,j)-1.d0)
               vfl(i,j)=1.d0-vfls(i,j)
            !   if(i==165.and.j==160) then
            !     print*,vfl(i,j),1.d0-vfls(i,j)+vofeps
            !     print*,'after test'
            !   end if
             end if
             if(vfl(i,j)<=1.d0-vfls(i,j)-vofeps.and.vfls(i,j)>epsi.and.        &
                vfl(i,j)/(1.d0-vfls(i,j))>=1.d0-vfl(i,j)*                      &
                dabs(dttol*((ue(i,j)-ue(i-1,j))/PGrid%dx(i,j)+                 &
                (ve(i,j)-ve(i,j-1))/PGrid%dy(i,j)))/(1.d0-vfls(i,j))) then
               sumweight=sumweight+(1.d0-vfl(i,j)-vfls(i,j))
             end if
           end do
         end do

         do j = 2,Jsize-1
           do i = 2,Isize-1
             if(vfl(i,j)<=1.d0-vofeps.and.vfls(i,j)<epsi.and.                  &
                vfl(i,j)>=1.d0-1.d-4.and.                                      &
               (vfls(i-1,j)>epsi.or.vfls(i+1,j)>epsi.or.                       &
                vfls(i,j-1)>epsi.or.vfls(i,j+1)>epsi)) then
               sumweight=sumweight+(1.d0-vfl(i,j))
             end if
           end do
         end do
         tempc=1
         ! If the residual of vfl in cells that excess 1.d0-vfls, is smaller
         ! than deficit of vfl in cells that need a 'treatment',
         ! more vfl from interface cell will contribute to this residual
         if(sumvol-sumweight<-vofeps) then
           voldif=sumweight-sumvol
           sumweight2=0.d0
           do j = 1,Jsize
             do i = 1,Isize
               if(vfl(i,j)<=0.6d0.and.vfls(i,j)<epsi.and.vfl(i,j)>0.4d0) then
                 sumweight2=sumweight2+(vfl(i,j)-0.5d0)
               end if
             end do
           end do
           do j = 1,Jsize
             do i = 1,Isize
               if(vfl(i,j)<=0.6d0.and.vfls(i,j)<epsi.and.vfl(i,j)>0.4d0) then
                 vfl(i,j)=vfl(i,j)-voldif*(vfl(i,j)-0.5d0)/sumweight2
               end if
             end do
           end do
           sumvol=sumvol+voldif
         end if
         ! Adding volume fraction of fluid to the cells that should contain only fluid
         ! However, numerical error create  small gas volume fraction inside this cell
         ! due to the VFLGas = 1.d0-vfl-vfls
         SumAllVolBig=0.d0
         SumAllVolSmall=0.d0
         do j = 2,Jsize-1
           do i = 2,Isize-1
             if(vfl(i,j)<=1.d0-vofeps.and.vfls(i,j)<epsi.and.                  &
               vfl(i,j)>=1.d0-1.d-4.and.                                       &
               (vfls(i-1,j)>epsi.or.vfls(i+1,j)>epsi.or.                       &
                vfls(i,j-1)>epsi.or.vfls(i,j+1)>epsi)) then
               voldif=sumvol*(1.d0-vfl(i,j))/sumweight
               vfl(i,j)=vfl(i,j)+voldif
             end if
           end do
         end do
         do j = 1,Jsize
           do i = 1,Isize
             if(vfl(i,j)<=1.d0-vfls(i,j)-vofeps.and.vfls(i,j)>epsi.and.        &
                vfl(i,j)/(1.d0-vfls(i,j))>=1.d0-vfl(i,j)*                      &
                dabs(dttol*((ue(i,j)-ue(i-1,j))/PGrid%dx(i,j)+                 &
               (ve(i,j)-ve(i,j-1))/PGrid%dy(i,j)))/(1.d0-vfls(i,j))) then
               voldif=sumvol*(1.d0-vfl(i,j)-vfls(i,j))/sumweight
               vfl(i,j)=vfl(i,j)+voldif
             end if
             if(vfl(i,j)>=1.d0-vfls(i,j)-vofeps) then
               SumAllVolBig=SumAllVolBig+(vfl(i,j)-1.d0+vofeps+vfls(i,j))
               vfl(i,j)=1.d0-vfls(i,j)
             end if
             if(vfl(i,j)<vofeps) then
               SumAllVolSmall=SumAllVolSmall+vofeps-vfl(i,j)
               vfl(i,j)=0.d0
             end if
           end do
         end do
         sumweight2=0.d0
         voldif=SumAllVolBig-SumAllVolSmall
         do j = 1,Jsize
           do i = 1,Isize
             if(vfl(i,j)<=1.d0-vfls(i,j)-1.d-5.and.vfls(i,j)<epsi.and.         &
                                                  vfl(i,j)>0.5d0) then
               sumweight2=sumweight2+(vfl(i,j)-0.5d0)
             end if
           end do
         end do
         do j = 1,Jsize
           do i = 1,Isize
             if(vfl(i,j)<=1.d0-vfls(i,j)-1.d-5.and.vfls(i,j)<epsi.and.         &
                                                  vfl(i,j)>0.5d0) then
               vfl(i,j)=vfl(i,j)+voldif*(vfl(i,j)-0.5d0)/sumweight2
             end if
             if(vfl(i,j)>=1.d0-vfls(i,j)-vofeps) vfl(i,j)=1.d0-vfls(i,j)
             if(vfl(i,j)<vofeps) vfl(i,j)=0.d0
           end do
         end do
         UpdateNorVec=1
         call Interface_Reconstruct(PGrid,nx,ny,dis,UpdateNorVec)
         PCell%nx(:,:) = nx(:,:)
         PCell%ny(:,:) = ny(:,:)
         call Redistance(PGrid,nx,ny,dis)
         do i=1,Isize
           do j=1,Jsize
             if(isnan(phi(i,j))) then
               print*,i,j
               print*,'after resistance'
             end if
           end do
         end do
    !     call Print_Result_Tecplot_VarP(PGrid,PCell%vof,INT8(k))
       end do
   !    print*,'test'
   !    print*,BoomCase%Posp%y,BoomCase%vs
   !    print*,BoomCase%asy,dt
   !    print*,BoomCase%Posp%y-BoomCase%asy*dt**2.d0/2.d0
   !    print*,
       nullify(vfl,phi,vflS,phiS,nxs,nys)
       deallocate(nx,ny)
       deallocate(dis)
       deallocate(temvfx)
       deallocate(temvfy)
       deallocate(temlsx)
       deallocate(temlsy)
       deallocate(ue)
       deallocate(ve)
    end subroutine Coupled_LS_VOF

    subroutine Boundary_Condition(vari)
      IMPLICIT NONE
      REAL(KIND=dp),DIMENSION(:,:)::vari
      INTEGER(kind=it4b):: i,j
      do i = 1,ISize
        vari(i,1) = vari(i,2)
        vari(i,JSize) = vari(i,JSize-1)
      end do
      do j = 1,JSize
        vari(1,j) = vari(2,j)
        vari(ISize,j) = vari(ISize-1,j)
      end do
    end subroutine Boundary_Condition

    subroutine Boundary_Condition_Vof(vari)
       IMPLICIT NONE
       REAL(KIND=dp),DIMENSION(:,:),INTENT(INOUT)::vari
       INTEGER:: i,j
       do i = 1,ISize
         vari(i,1)=dmin1(vari(i,2)+vfls(i,1),1.d0)-vfls(i,1)
         vari(i,JSize)=dmin1(vari(i,JSize-1)+vfls(i,Jsize),1.d0)-vfls(i,Jsize)
       end do
       do j = 1,JSize
         vari(1,j)=dmin1(vari(2,j)+vfls(1,j),1.d0)-vfls(1,j)
         vari(ISize,j)=dmin1(vari(Isize-1,j)+vfls(Isize,j),1.d0)-vfls(ISize,j)
       end do
    end subroutine Boundary_Condition_Vof

    SUBROUTINE X_Sweep(PGrid,PCell,temvf,temls,ue,ve,nx,ny,dis,dtv,itt,        &
                                         UpdateNorVec,NondiT,k)
       IMPLICIT NONE
       TYPE(Grid)                                            :: PGrid
       TYPE(Cell)                                            :: PCell
       INTEGER(kind=it4b),INTENT(IN)                         :: UpdateNorVec,k
       INTEGER(kind=it8b),INTENT(IN)                         :: itt
       REAL(KIND=dp),INTENT(IN)                              :: dtv,NondiT
       REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(IN)   :: ue,ve
       REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(INOUT):: nx,ny,dis
       REAL(KIND=dp),DIMENSION(:,:),allocatable              :: temvf,temls
       INTEGER(kind=it4b)                                    :: i,j
       REAL(KIND=dp)                                         :: flux,lsr,etau,ul1,xwu
       call Interface_Reconstruct(PGrid,nx,ny,dis,UpdateNorVec)
       flux = 0.d0
     ! volume of fluid
       xwu=PGrid%x(1,1)-PGrid%dx(1,1)/2.d0
    !   if(itt==361) print*,vfl(1,155)
       do j = 1,JSize
         do i = 1,ISize
           if(vflS(i,j)<1.d0-vofeps) temvf(i,j)=temvf(i,j)/(1.d0-vflS(i,j))
         end do
       end do
       do j = 1,JSize
         do i = 1,ISize
           if(ue(i,j)>=0.d0) then
             if(vfl(i,j)>=(1.d0-vofeps).or.vfl(i,j)<=vofeps)then
               flux = vfl(i,j)*ue(i,j)*dtv
               if(vfl(i,j)>=1.d0) flux=ue(i,j)*dtv
               if(vfl(i,j)<0.d0) flux=0.d0
             else
               call rightflux(nx(i,j),ny(i,j),dis(i,j),vfl(i,j),nxs(i,j),      &
                    nys(i,j),phiS(i,j),vfls(i,j),PGrid%dx(i,j),PGrid%dy(i,j),  &
                                                            ue(i,j)*dtv,flux)
             end if
             if(PCell%EEdge_Area(i,j)>1.d0-vofeps.and.                         &
                                 vfls(i,j)>1.d0-vofeps) flux=1.d0*ue(i,j)*dtv
           else
             if(i<Isize) then
               if(vfl(i+1,j)>=(1.d0-vofeps).or.vfl(i+1,j)<=vofeps) then
                 flux=vfl(i+1,j)*ue(i,j)*dtv
                 if(vfl(i+1,j)>=1.d0) flux=ue(i,j)*dtv
                 if(vfl(i+1,j)<0.d0) flux=0.d0
               else
                 call leftflux(nx(i+1,j),ny(i+1,j),dis(i+1,j),vfl(i+1,j),      &
                           nxs(i+1,j),nys(i+1,j),phis(i+1,j),vfls(i+1,j),      &
                           PGrid%dx(i+1,j),PGrid%dy(i+1,j),-ue(i,j)*dtv,flux)
                 flux=-flux
               end if
               if(PCell%WEdge_Area(i+1,j)>1.d0-vofeps.and.                     &
                                 vfls(i+1,j)>1.d0-vofeps) flux=1.d0*ue(i,j)*dtv
             else
               Flux=vfl(i,j)*ue(i,j)*dtv
             end if
           end if
           if(i>=1)temvf(i,j)=temvf(i,j)-flux*PCell%EEdge_Area(i,j)/PGrid%dx(i,j)
           if(i<=ISize-1) then
             temvf(i+1,j)=temvf(i+1,j)+flux*PCell%WEdge_Area(i+1,j)/PGrid%dx(i+1,j)
           end if
           if(isnan(temvf(i,j))) then
             print*,flux
             print*,'X-Sweep 1263'
           end if
         end do
     !   For boundary condition
         etau=amp0*dsin(kw*xwu-omew*(NondiT+dble(k)*dtv))
     !    UwInlet-Amp0*kw*(UwInlet-cw0)*                   &
     !               dsin(kw*(xwu-cw0*Time))*dcosh(kw*PGrid%y(1,j))/dsinh(kw*Hw)
         if(PGrid%y(1,j)-PGrid%dy(1,j)/2.d0<Hw+etau) then
           ul1=ue(0,j)!kw*cw0*etau*dcosh(kw*PGrid%y(1,j))/dsinh(kw*Hw)!
     !      Flux=dmin1(1.d0,dabs(Hw+etau-(PGrid%y(1,j)-PGrid%dy(1,j)/2.d0))/    &
     !                  PGrid%dy(1,j))*ul1*dtv
           Flux=vfl(1,j)*ul1*dtv
         else
           Flux=0.d0
         end if
         temvf(1,j)=temvf(1,j)+flux*PCell%WEdge_Area(1,j)/PGrid%dx(1,j)
       end do
     ! level set
       lsr = 0.d0
       flux = 0.d0
       do j = 1,JSize
         do i = 2,ISize-1
           if(ue(i,j)>=0.d0) then
          !   lsr =  phi(i,j)+del_x/2.d0*(1.d0-ue(i,j)*dtv/del_x)*              &
          !         (phi(i+1,j)-phi(i-1,j))/(2.d0*del_x)
              lsr=(phi(i+1,j)+phi(i-1,j))/2.d0
           else
             if(i<=ISize-2) then
           !    lsr = phi(i+1,j)-del_x/2.d0*(1.d0+ue(i,j)*dtv/del_x)*           &
           !          (phi(i+2,j)-phi(i,j))/(2.d0*del_x)
               lsr=(phi(i+2,j)+phi(i,j))/2.d0
             else
           !    lsr = phi(i+1,j)-del_x/2.d0*(1.d0+ue(i,j)*dtv/del_x)*           &
           !          (phi(i+1,j)-phi(i,j))/(del_x)
               lsr=(phi(i+1,j)+phi(i-1,j))/2.d0
             end if
           end if
           flux=ue(i,j)*lsr*dtv
           if(i>=2) temls(i,j)=temls(i,j)-flux*PCell%EEdge_Area(i,j)/PGrid%dx(i,j)
           if(i<=ISize-1) then
             temls(i+1,j)=temls(i+1,j)+flux*PCell%WEdge_Area(i+1,j)/PGrid%dx(i+1,j)
           end if
         end do
       end do
       do j=1,Jsize
         if(ue(1,j)>=0.d0) then
           lsr=phi(1,j)+PGrid%dx(1,j)/2.d0*(1.d0-ue(1,j)*dtv/PGrid%dx(1,j))*   &
               (phi(2,j)-phi(1,j))/(PGrid%x(2,j)-PGrid%x(1,j))
         else
           lsr=phi(2,j)-PGrid%dx(2,j)/2.d0*(1.d0+ue(1,j)*dtv/PGrid%dx(2,j))*   &
               (phi(3,j)-phi(1,j))/(2.d0*PGrid%dx(2,j))
         end if
         flux=lsr*ue(1,j)*dtv
         temls(1,j)=temls(1,j)-flux*PCell%EEdge_Area(1,j)/PGrid%dx(1,j)
         temls(2,j)=temls(2,j)+flux*PCell%EEdge_Area(1,j)/PGrid%dx(2,j)
         temls(1,j)=temls(1,j)+phi(1,j)*ue(0,j)*dtv*PCell%WEdge_Area(1,j)/     &
                                                          PGrid%dx(1,j)
         if(ue(Isize,j)>0.d0) then
           lsr=phi(Isize,j)+PGrid%dx(Isize,j)/2.d0*(1.d0-ue(Isize,j)*dtv/      &
               PGrid%dx(Isize,j))*(phi(Isize,j)-phi(Isize-1,j))/               &
               (PGrid%x(Isize,j)-PGrid%x(Isize-1,j))
         else
           lsr=phi(Isize,j)-PGrid%dx(Isize-1,j)/2.d0*(1.d0+ue(Isize,j)*dtv/    &
               PGrid%dx(Isize-1,j))*(phi(Isize,j)-phi(Isize-1,j))/             &
               (PGrid%x(Isize,j)-PGrid%x(Isize-1,j))
         end if
         flux=lsr*ue(Isize,j)*dtv
         temls(Isize,j)=temls(Isize,j)-flux*PCell%EEdge_Area(Isize,j)/         &
                                                            PGrid%dx(Isize,j)
       end do
    end subroutine X_Sweep

    subroutine Y_Sweep(PGrid,PCell,temvf,temls,ue,ve,nx,ny,dis,dtv,itt,        &
                                                                 UpdateNorVec)
       IMPLICIT NONE
       TYPE(Grid),INTENT(IN):: PGrid
       TYPE(Cell),INTENT(IN):: PCell
       INTEGER(kind=it4b),INTENT(IN):: UpdateNorVec
       INTEGER(kind=it8b),INTENT(IN):: itt
       REAL(KIND=dp),INTENT(IN):: dtv
       REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(IN):: ue,ve
       REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(INOUT):: nx,ny,dis
       REAL(KIND=dp),DIMENSION(:,:),allocatable:: temvf,temls
       INTEGER(kind=it4b):: i,j
       REAL(KIND=dp):: flux,lst
       call Interface_Reconstruct(PGrid,nx,ny,dis,UpdateNorVec)
       flux = 0.d0
       ! volume of fluid
       do j = 1,JSize
         do i = 1,ISize
           if(vflS(i,j)<1.d0-vofeps) temvf(i,j)=temvf(i,j)/(1.d0-vflS(i,j))
         end do
       end do
       do i = 1,ISize
         do j = 1,JSize
           if(ve(i,j)>0.d0) then
             if(vfl(i,j)>=(1.d0-vofeps).or.vfl(i,j)<=vofeps) then
               flux=vfl(i,j)*ve(i,j)*dtv
               if(vfl(i,j)>=1.d0) flux=ve(i,j)*dtv
               if(vfl(i,j)<0.d0) flux=0.d0
             else
               call topflux(nx(i,j),ny(i,j),dis(i,j),vfl(i,j),nxs(i,j),        &
                    nys(i,j),phis(i,j),vfls(i,j),PGrid%dx(i,j),PGrid%dy(i,j),  &
                                                           ve(i,j)*dtv,flux)
             end if
             if(PCell%NEdge_Area(i,j)>1.d0-vofeps.and.                   &
                                 vfls(i,j)>1.d0-vofeps) flux=1.d0*ve(i,j)*dtv
           else
             if(j<Jsize) then
               if(vfl(i,j+1)>=(1.d0-vofeps).or.vfl(i,j+1)<=vofeps) then
                 flux=vfl(i,j+1)*ve(i,j)*dtv
                 if(vfl(i,j+1)>=1.d0) flux=ve(i,j)*dtv
                 if(vfl(i,j+1)<0.d0) flux=0.d0
               else
                 call bottomflux(nx(i,j+1),ny(i,j+1),dis(i,j+1),vfl(i,j+1),    &
                             nxs(i,j+1),nys(i,j+1),phis(i,j+1),vfls(i,j+1),    &
                             PGrid%dx(i,j+1),PGrid%dy(i,j+1),-ve(i,j)*dtv,flux)
                 flux=-flux
               endif
               if(PCell%SEdge_Area(i,j+1)>1.d0-vofeps.and.                     &
                                 vfls(i,j+1)>1.d0-vofeps) flux=1.d0*ve(i,j)*dtv
             else
                flux=vfl(i,j)*ve(i,j)*dtv
             endif
           end if
           if(j>=1)temvf(i,j)=temvf(i,j)-flux*PCell%NEdge_Area(i,j)/PGrid%dy(i,j)
           if(j<=JSize-1) then
             temvf(i,j+1)=temvf(i,j+1)+flux*PCell%SEdge_Area(i,j+1)/PGrid%dy(i,j+1)
           end if
        !   if(itt==45.and.i==299.and.(j==124.or.j==123)) then
        !     print*,'Inside Y-Sweep'
        !     print*,flux
        !     print*,ve(i,j)*dtv
        !     print*,vfl(i,j+1)+vfls(i,j+1)
        !     print*,
        !   end if
         end do
         temvf(i,1)=temvf(i,1)+VofInlet*ve(i,0)*dtv*PCell%SEdge_Area(i,1)/     &
                                                                  PGrid%dy(i,1)
       end do

       lst = 0.d0
       flux = 0.d0
    ! level set
       do j = 2,JSize-1
         do i = 1,ISize
           if(ve(i,j)>=0.d0) then
        !     lst = phi(i,j)+del_y/2.d0*(1.d0-ve(i,j)*dtv/del_y)*               &
        !           (phi(i,j+1)-phi(i,j-1))/(2.d0*del_y)
             lst=(phi(i,j+1)+phi(i,j-1))/2.d0
           else
             if(j<=JSize-2) then
          !     lst = phi(i,j+1)-del_y/2.d0*(1.d0+ve(i,j)*dtv/del_y)*           &
          !          (phi(i,j+2)-phi(i,j))/(2.d0*del_y)
               lst=(phi(i,j+2)+phi(i,j))/2.d0
             else
           !    lst = phi(i,j+1)-del_y/2.d0*(1.d0+ve(i,j)*dtv/del_y)*           &
           !    (phi(i,j+1)-phi(i,j))/del_y
               lst=(phi(i,j+1)+phi(i,j-1))/2.d0
             end if
           end if
           flux = lst*ve(i,j)*dtv
           if(j>=2) temls(i,j)=temls(i,j)-flux*PCell%NEdge_Area(i,j)/PGrid%dy(i,j)
           if(j<=JSize-1) then
             temls(i,j+1)=temls(i,j+1)+flux*PCell%SEdge_Area(i,j+1)/PGrid%dy(i,j+1)
           end if
         end do
       end do

       do i=1,Isize
         if(ve(i,1)>=0.d0) then
           lst=phi(i,1)+PGrid%dy(i,1)/2.d0*(1.d0-ve(i,1)*dtv/PGrid%dy(i,1))*   &
               (phi(i,2)-phi(i,1))/(PGrid%y(i,2)-PGrid%y(i,1))
         else
           lst=phi(i,2)-PGrid%dy(i,2)/2.d0*(1.d0+ve(i,1)*dtv/PGrid%dy(i,2))*   &
                    (phi(i,3)-phi(i,1))/(PGrid%y(i,3)-PGrid%y(i,1))
         end if
         flux = lst*ve(i,1)*dtv
         temls(i,1)=temls(i,1)-flux*PCell%NEdge_Area(i,1)/PGrid%dy(i,1)
         temls(i,2)=temls(i,2)+flux*PCell%NEdge_Area(i,1)/PGrid%dy(i,2)
         temls(i,1)=temls(i,1)+phi(i,1)*ve(i,0)*dtv*PCell%SEdge_Area(i,1)/     &
                                                                PGrid%dy(i,1)
         if(ve(i,Jsize)>0.d0) then
           lst=phi(i,Jsize)+PGrid%dy(i,Jsize)/2.d0*(1.d0-ve(i,Jsize)*dtv/      &
               PGrid%dy(i,Jsize))*(phi(i,Jsize)-phi(i,Jsize-1))/               &
               (PGrid%y(i,Jsize)-PGrid%y(i,Jsize-1))
         else
           lst=phi(i,Jsize)-PGrid%dy(i,Jsize-1)/2.d0*(1.d0+ve(i,Jsize)*dtv/    &
               PGrid%dy(i,Jsize-1))*(phi(i,Jsize)-phi(i,Jsize-1))/             &
               (PGrid%y(i,Jsize)-PGrid%y(i,Jsize-1))
         end if
         flux=lst*ve(i,Jsize)*dtv
         temls(i,Jsize)=temls(i,Jsize)-flux/PGrid%dy(i,Jsize)
       end do
    end subroutine Y_Sweep

  ! level set step
  ! build-up interface

    subroutine Interface_Reconstruct(PGrid,nx,ny,dis,UpdateNorVec)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid
      INTEGER(kind=it4b),INTENT(IN):: UpdateNorVec
      REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(INOUT)::nx,ny,dis
      REAL(KIND=dp):: nxx,nyy,diss,temp
      INTEGER:: i,j,flag
      do i = 1,ISize
        do j = 1,JSize
          flag=0
          if(vfl(i,j)>=vofeps.and.vfl(i,j)<(1.d0-vofeps-vfls(i,j)).and.            &
                                                        vfls(i,j)<epsi)then
            flag=1
          end if
          if(vfl(i,j)>=vofeps.and.vfl(i,j)<(1.d0-vofeps-vfls(i,j)).and.            &
                                vfls(i,j)>=epsi.and.vfls(i,j)<1.d0-epsi)then
            flag=2
          end if
          if(flag==1.or.flag==2) then
            if(UpdateNorVec==1) then
              if(i>1.and.i<Isize.and.j>1.and.j<Jsize) then
                call Normal_Vector_Irre(PGrid,i,j,nxx,nyy)
              else
                if(i==1)then
                  nxx=(phi(i+1,j)-phi(i,j))/(PGrid%x(i+1,j)-PGrid%x(i,j))
                elseif(i==Isize)then
                  nxx=(phi(i,j)-phi(i-1,j))/(PGrid%x(i,j)-PGrid%x(i-1,j))
                endif
                if(j==1)then
                  nyy=(phi(i,j+1)-phi(i,j))/(PGrid%y(i,j+1)-PGrid%y(i,j))
                elseif(j==Jsize)then
                  nyy=(phi(i,j)-phi(i,j-1))/(PGrid%y(i,j)-PGrid%y(i,j-1))
                endif
              endif
          !   Test the normal vector for cell containing all three phases
              if(vfls(i,j)>epsi.and.vfls(i,j)<1.d0-epsi) then
          !      nxx=dsign(1.d0,nxs(i,j))*dabs(nxx)
              end if
              if(vfl(i,j)>epsi.and.vfl(i,j)<1.d0-epsi) then
          !      nxx=dsign(1.d0,nxs(i,j))*dabs(nxx)
              end if
          !   *******************************************************************
              temp = dsqrt(nxx**2.d0+nyy**2.d0)
              if(temp<1.d-14) then
                nxx = 0.d0
                nyy = 1.d0
              else
                nxx = nxx/temp
                nyy = nyy/temp
              end if
            end if
            select case(flag)
              case(1) ! for only 2 phases (air(gas) and liquid(fluid))
                call distance(i,j,PGrid%dx(i,j),PGrid%dy(i,j),nxx,nyy,diss)
              case(2)
                call DistanceFluidCalculate(nxs(i,j),nys(i,j),phis(i,j),nxx,   &
                                nyy,vfl(i,j),PGrid%dx(i,j),PGrid%dy(i,j),diss)
            end select
          else
            nxx=0.d0;nyy=1.d0
            if(1.d0-vfls(i,j)>vofeps) then
              diss=(0.5d0-vfl(i,j)/(1.d0-vfls(i,j)))*PGrid%dy(i,j)
            else
              diss=(0.d0-vfl(i,j))*PGrid%dy(i,j)
            end if
          end if
          nx(i,j) = nxx
          ny(i,j) = nyy
          dis(i,j) = diss
        end do
      end do
    end subroutine Interface_Reconstruct

    subroutine distance(i,j,dx,dy,nxx,nyy,diss)
       IMPLICIT NONE
       INTEGER i,j
       real(dp) slop_eps,nxx,nyy,nx,ny,diss,equal_eps
       real(dp) dll,dul,dlr,dur,vful,vflr,theta,dx,dy
       slop_eps = 1.d-20
       equal_eps = 1.d-20
       nx = dabs(nxx)
       ny = dabs(nyy)
       if(nx < slop_eps) then
          diss = (0.5d0-vfl(i,j))*dy
          return
       end if
       if(ny < slop_eps) then
          diss = (0.5d0-vfl(i,j))*dx
          return
       end if
       dll = 0.5d0*(nx*dx+ny*dy)
       dul = 0.5d0*(nx*dx-ny*dy)
       dlr = -dul
       dur = -dll
       call frac(nx,ny,dul,dx,dy,vful)
       call frac(nx,ny,dlr,dx,dy,vflr)
       if(vfl(i,j)<=vful.and.vfl(i,j)<=vflr) then
          if(vful>vflr) then
             call swap(vful,vflr)
             call swap(dul,dlr)
          end if
          diss = dll-dsqrt(((dul-dll)**2.d0)*vfl(i,j)/vful)
          return
       elseif(vfl(i,j)>=vful.and.vfl(i,j)>=vflr) then
          if(vflr<vful) then
             call swap(vful,vflr)
             call swap(dul,dlr)
          end if
          diss = dur+dsqrt((dlr-dur)**2.d0*(1.d0-vfl(i,j))/(1.d0-vflr))
          return
       else
          if(dabs(vflr-vful)<equal_eps)THEN
	        diss=0.5d0*(dul+dlr)
	        return
	      end if
	      theta=(vfl(i,j)-vful)/(vflr-vful)
	      diss=dlr*theta+dul*(1.0d0-theta)
	   end if
	   return
	end subroutine

    subroutine frac(nx,ny,diss,dx,dy,vrt)
       IMPLICIT NONE
       real(dp) slop_eps,nx,ny,diss,vrt,tnx,tny,dx,dy
       real(dp) xx,yy,topvf,rightvf,totalarea
       slop_eps = 1.d-20
       tnx = dabs(nx)
       tny = dabs(ny)
       if(diss+0.5d0*(tnx*dx+tny*dy)<=0.d0) then
          vrt = 1.d0
          return
       end if
       if(diss-0.5d0*(tnx*dx+tny*dy)>=0.d0) then
          vrt = 0.d0
          return
       end if
       if(tnx <= slop_eps) then
          vrt = 0.5d0-diss/dy
          return
       end if
       if(tny <= slop_eps) then
          vrt = 0.5d0-diss/dx
          return
       end if
       xx = (0.5d0*tny*dy-diss)/dx/tnx+0.5d0
       yy = (0.5d0*tnx*dx-diss)/dy/tny+0.5d0
       totalarea = 0.5d0*xx*yy
       topvf = 0.d0
       if(yy > 1.d0) topvf =((yy-1.d0)/yy)**2.d0
       rightvf = 0.d0
       if(xx>1.d0) rightvf = ((xx-1.d0)/xx)**2.d0
       vrt=totalarea*(1.d0-topvf-rightvf)
       if(vrt<0.d0) write(*,*) 'vrt less than 0'
       return
    end subroutine

    subroutine swap(a,b)
       IMPLICIT NONE
       real(dp) a,b,t
       t = a
       a = b
       b = t
       return
    end subroutine

!   calculate flux through the top of cell
    subroutine topflux(nxx,nyy,diss,volf,nxss,nyss,phiss,vols,dx,dy,vdt,flux)
       IMPLICIT NONE
       REAL(KIND=dp),INTENT(IN):: nxx,nyy,diss,volf,nxss,nyss,phiss,vols,vdt,dx,dy
       REAL(KIND=dp),INTENT(OUT):: flux
       REAL(KIND=dp):: eps,volss
       eps = 1.d-14
       if(vdt==0.d0) then
          flux = 0.d0
          return
       end if
       if(vols<epsi) then
         call frac(nxx,nyy,diss-0.5d0*(vdt-dy)*nyy,dx,vdt,flux)
       elseif(vols>=1.d0-epsi) then
         flux=volf
         return
       else
         call CellGeoCal(nxss,nyss,phiss-0.5d0*(vdt-dy)*nyss,nxx,nyy,          &
                                    diss-0.5d0*(vdt-dy)*nyy,dx,vdt,volss,flux)
         if(volf/(1.d0-vols)>=1.d0-vofeps) flux=1.d0-volss
       ! Using the volume averaging value
         flux=flux/(1.d0-volss+toldeno)
         if(volf+vols>1.d0-eps) flux=1.d0
       end if
       flux=flux*vdt
       return
    end subroutine

!   calculate flux through the bottom of cell
    subroutine bottomflux(nxx,nyy,diss,volf,nxss,nyss,phiss,vols,dx,dy,vdt,flux)
       IMPLICIT NONE
       REAL(KIND=dp),INTENT(IN):: nxx,nyy,diss,volf,nxss,nyss,phiss,vols,vdt,dx,dy
       REAL(KIND=dp),INTENT(OUT):: flux
       REAL(KIND=dp):: eps,volss
       eps = 1.d-14
       if(vdt==0.d0) then
          flux = 0.d0
          return
       end if
       if(vols<epsi) then
         call frac(nxx,nyy,diss+0.5d0*(vdt-dy)*nyy,dx,vdt,flux)
       elseif(vols>=1.d0-epsi) then
         flux=volf
         return
       else
         call CellGeoCal(nxss,nyss,phiss+0.5d0*(vdt-dy)*nyss,nxx,nyy,          &
                                    diss+0.5d0*(vdt-dy)*nyy,dx,vdt,volss,flux)
         if(volf/(1.d0-vols)>=1.d0-vofeps)  flux=1.d0-volss
       ! Using the volume averaging value
         flux=flux/(1.d0-volss+toldeno)
         if(volf+vols>1.d0-eps) flux=1.d0
       end if
       flux=flux*vdt
    end subroutine

!   calculate flux through the right of cell
    subroutine rightflux(nxx,nyy,diss,volf,nxss,nyss,phiss,vols,dx,dy,udt,flux)
       IMPLICIT NONE
       REAL(KIND=dp),INTENT(IN):: nxx,nyy,diss,volf,nxss,nyss,phiss,vols,udt,dx,dy
       REAL(KIND=dp),INTENT(OUT):: flux
       REAL(KIND=dp):: eps,volss
       eps = 1.d-14
       if(udt==0.d0) then
         flux = 0.d0
         return
       end if
       if(vols<epsi) then
         call frac(nxx,nyy,diss-0.5d0*(udt-dx)*nxx,udt,dy,flux)
       elseif(vols>=1.d0-epsi) then
         flux=volf
         return
       else
         call CellGeoCal(nxss,nyss,phiss-0.5d0*(udt-dx)*nxss,nxx,nyy,          &
                                    diss-0.5d0*(udt-dx)*nxx,udt,dy,volss,flux)
         if(volf/(1.d0-vols)>=1.d0-vofeps)  flux=1.d0-volss
       ! Using the volume averaging value
         flux=flux/(1.d0-volss+toldeno)
         if(volf+vols>1.d0-eps) flux=1.d0
       end if
       flux = flux*udt
    end subroutine

! calculate flux throught the left of cell
    subroutine leftflux(nxx,nyy,diss,volf,nxss,nyss,phiss,vols,dx,dy,udt,flux)
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(IN):: nxx,nyy,diss,volf,nxss,nyss,phiss,vols,udt,dx,dy
      REAL(KIND=dp),INTENT(OUT):: flux
      REAL(KIND=dp):: eps,volss
      eps = 1.d-14
      if(udt==0.d0) then
        flux = 0.d0
        return
      end if
      if(vols<epsi) then
        call frac(nxx,nyy,diss+0.5d0*(udt-dx)*nxx,udt,dy,flux)
      elseif(vols>=1.d0-epsi) then
        flux=volf
        return
      else
        call CellGeoCal(nxss,nyss,phiss+0.5d0*(udt-dx)*nxss,nxx,nyy,           &
                                   diss+0.5d0*(udt-dx)*nxx,udt,dy,volss,flux)
        if(volf/(1.d0-vols)>=1.d0-vofeps) flux=1.d0-volss
      ! Using the volume averaging value
        flux=flux/(1.d0-volss+toldeno)
        if(volf+vols>1.d0-eps) flux=1.d0
      end if
      flux = flux*udt
    end subroutine

  ! Redistance from vof to level set
    subroutine Redistance(PGrid,nxx1,nyy1,diss1)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid
      REAL(KIND=dp),DIMENSION(:,:),allocatable,INTENT(IN):: nxx1,nyy1,diss1
      INTEGER(kind=it4b):: i,j,ii,jj,l,m
      REAL(KIND=dp):: xv,yv,dv,deuc,dij1,dvs,del_x,del_y
      REAL(KIND=dp):: xp,yp,xoff,yoff,xfc,yfc,xs,ys,tol
      logical:: pointp
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: phiaux
      INTEGER,DIMENSION(:,:),allocatable:: fix
      allocate(fix(Isize,Jsize))
      allocate(phiaux(Isize,Jsize))
      fix(:,:)=0
      phiaux(:,:)=1.d4
      tol=1.d-24
      del_x=PGrid%dx(Isize/2,Jsize/2)
      del_y=PGrid%dy(Isize/2,Jsize/2)
      do i = 1,ISize
        do j = 1,JSize
          if((0.d0+vofeps)<vfl(i,j).and.vfl(i,j)<(1.d0-vofeps-vflS(i,j))) then
            phiaux(i,j)=diss1(i,j)
            fix(i,j)=1
            do ii=-band_width,band_width
              do jj=-band_width,band_width
              ! determine the point xv on the boundary cell(i,j,k) with the
              ! shortest distance to the cell center of (i+-ii,j+-jj)
                if(i+ii>=1.and.i+ii<=ISize.and.j+jj>=1.and.j+jj<=JSize) then
                  if(vfl(i+ii,j+jj)<0.d0+vofeps.or.                              &
                    vfl(i+ii,j+jj)>=1.d0-vofeps-vflS(i+ii,j+jj)) then
                    fix(i+ii,j+jj) = 1
                    l=max(-1,min(1,ii))
                    m=max(-1,min(1,jj))
                    xv=PGrid%dx(i,j)*dble(l)/2.d0
                    yv=PGrid%dy(i,j)*dble(m)/2.d0
                    dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                    if(dvs<0.d0) then
                      if(dabs(yv)<vofeps) then
                        yv=PGrid%dy(i,j)/2.d0
                        dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                        if(dvs<0.d0) then
                          yv=-PGrid%dy(i,j)/2.d0
                          dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                          if(dvs<0.d0) then
                            xv=PGrid%dx(i,j)*(0.d0-dble(l))/2.d0
                            dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                            if(dvs<0.d0) yv=PGrid%dy(i,j)/2.d0
                          end if
                        end if
                      elseif(dabs(xv)<vofeps) then
                        xv=PGrid%dx(i,j)/2.d0
                        dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                        if(dvs<0.d0) then
                          xv=-PGrid%dx(i,j)/2.d0
                          dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                          if(dvs<0.d0) then
                            yv=PGrid%dy(i,j)*(0.d0-dble(m))/2.d0
                            dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                            if(dvs<0.d0) xv=PGrid%dx(i,j)/2.d0
                          end if
                        end if
                      else
                        yv=PGrid%dy(i,j)*(0.d0-dble(m))/2.d0
                        dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                        if(dvs<0.d0) then
                          xv=PGrid%dx(i,j)*(0.d0-dble(l))/2.d0
                          yv=PGrid%dy(i,j)*dble(m)/2.d0
                          dvs=nxs(i,j)*xv+nys(i,j)*yv+phiS(i,j)
                          if(dvs<0.d0) yv=PGrid%dy(i,j)*(0.d0-dble(m))/2.d0
                        end if
                      end if
                    end if

                    dv=nxx1(i,j)*xv+nyy1(i,j)*yv+diss1(i,j)
                    if(dv*dsign(1.d0,phi(i+ii,j+jj))<=0.d0) then
                      deuc=dsqrt(((PGrid%x(i+ii,j+jj)-PGrid%x(i,j))-xv)**2.d0+ &
                                 ((PGrid%y(i+ii,j+jj)-PGrid%y(i,j))-yv)**2.d0)
                      if(vfls(i+ii,j+jj)<epsi) then
                        phiaux(i+ii,j+jj)=dsign(1.d0,0.5d0-vfl(i+ii,j+jj))*    &
                                          dmin1(deuc,dabs(phiaux(i+ii,j+jj)))
                      else
                        if(dabs(deuc)<dabs(phiaux(i+ii,j+jj)).and.             &
                                  vfl(i+ii,j+jj)>=1.d0-vofeps-vflS(i+ii,j+jj)) &
                          phiaux(i+ii,j+jj)=-dabs(deuc)
                        if(dabs(deuc)<dabs(phiaux(i+ii,j+jj)).and.             &
                                            vfl(i+ii,j+jj)<vofeps)             &
                          phiaux(i+ii,j+jj)=dabs(deuc)
                      endif
                   ! third step: find the projection of x' onto the interface
                    else
                      dij1=nxx1(i,j)*(PGrid%x(i+ii,j+jj)-PGrid%x(i,j))+        &
                          nyy1(i,j)*(PGrid%y(i+ii,j+jj)-PGrid%y(i,j))+diss1(i,j)
                      xp=(PGrid%x(i+ii,j+jj)-PGrid%x(i,j))-dij1*nxx1(i,j)
                      yp=(PGrid%y(i+ii,j+jj)-PGrid%y(i,j))-dij1*nyy1(i,j)
                      pointp=isinsiderect(xp,yp,PGrid%dx(i,j),PGrid%dy(i,j))
                      if((pointp.eqv..true.).and.vfls(i+ii,j+jj)<epsi) then
                        deuc=dsqrt(((PGrid%x(i+ii,j+jj)-PGrid%x(i,j))-xp)**2.d0+&
                                   ((PGrid%y(i+ii,j+jj)-PGrid%y(i,j))-yp)**2.d0) ! deuc = dij1
                        phiaux(i+ii,j+jj)=dsign(1.d0,0.5d0-vfl(i+ii,j+jj))*    &
                                          dmin1(deuc,dabs(phiaux(i+ii,j+jj)))
                      elseif((pointp.eqv..true.).and.vfls(i+ii,j+jj)>=epsi.and. &
                                            vfls(i+ii,j+jj)<1.d0-epsi) then
                        deuc=dsqrt(((PGrid%x(i+ii,j+jj)-PGrid%x(i,j))-xp)**2.d0+&
                                   ((PGrid%y(i+ii,j+jj)-PGrid%y(i,j))-yp)**2.d0)
                        if(dabs(deuc)<dabs(phiaux(i+ii,j+jj)).and.             &
                                   vfl(i+ii,j+jj)>=1.d0-vofeps-vfls(i+ii,j+jj))&
                          phiaux(i+ii,j+jj)=-dabs(deuc)
                        if(dabs(deuc)<dabs(phiaux(i+ii,j+jj)).and.             &
                                            vfl(i+ii,j+jj)<vofeps)             &
                          phiaux(i+ii,j+jj)=dabs(deuc)
                      elseif(vfls(i+ii,j+jj)>=1.d0-epsi) then
                        if(vfls(i,j)<epsi) then
                          deuc=dsqrt(((PGrid%x(i+ii,j+jj)-PGrid%x(i,j))-xp)**2.d0+&
                                     ((PGrid%y(i+ii,j+jj)-PGrid%y(i,j))-yp)**2.d0)
                          if(dabs(deuc)<dabs(phiaux(i+ii,j+jj)))               &
                             phiaux(i+ii,j+jj)=dsign(1.d0,dij1)*deuc
                        end if
                      else
                    ! forth step: find the corner of
                        xoff=dmax1(dabs(xp)-0.5d0*PGrid%dx(i,j),0.d0)
                        yoff=dmax1(dabs(yp)-0.5d0*PGrid%dy(i,j),0.d0)
                        xfc=dsign(1.d0,xp)*0.5d0*PGrid%dx(i,j)
                        yfc=dsign(1.d0,yp)*0.5d0*PGrid%dy(i,j)
                        if(xoff*dabs(nxx1(i,j)+tol)>=yoff*dabs(nyy1(i,j)+tol))then
                          xs=xfc
                          ys=(diss1(i,j)+nxx1(i,j)*xs)/(-nyy1(i,j)+tol)
                        else
                          ys=yfc
                          xs=(diss1(i,j)+nyy1(i,j)*ys)/(-nxx1(i,j)+tol)
                        endif
                        deuc=dsqrt(((PGrid%x(i+ii,j+jj)-PGrid%x(i,j))-xs)**2.d0+&
                                   ((PGrid%y(i+ii,j+jj)-PGrid%y(i,j))-ys)**2.d0)
                        dvs=nxs(i,j)*xs+nys(i,j)*ys+phis(i,j)
                        if(dvs>=0.d0) then
                          if(vfls(i+ii,j+jj)<epsi) then
                            phiaux(i+ii,j+jj)=dsign(1.d0,0.5d0-vfl(i+ii,j+jj))*&
                                            dmin1(deuc,dabs(phiaux(i+ii,j+jj)))
                          else
                            if(dabs(deuc)<dabs(phiaux(i+ii,j+jj)).and.         &
                               vfl(i+ii,j+jj)>=1.d0-vofeps-vfls(i+ii,j+jj))    &
                            phiaux(i+ii,j+jj)=-dabs(deuc)
                            if(dabs(deuc)<dabs(phiaux(i+ii,j+jj)).and.         &
                                                 vfl(i+ii,j+jj)<vofeps)        &
                            phiaux(i+ii,j+jj)=dabs(deuc)
                          end if
                        end if
                      endif
                    endif
                  endif
                  ! this condition will eliminate bugs in the case that
                  ! interface is vertical or horizontal
                  if(phiaux(i+ii,j+jj)>9.999d3) fix(i+ii,j+jj)=0
                end if
              end do
            end do
          end if
        end do
      end do
      do i = 1,ISize
        do j = 1,JSize
          if(fix(i,j)==1) then
            if(phiaux(i,j)>1.d3) then
              print*,i,j,phi(i,j),vfl(i,j),vfls(i,j)
              print*,phiaux(i,j)
              pause 'fuck you bugs, Clsvof_Mod 1053'
            end if
            phi(i,j)=phiaux(i,j)
          else
            if(vfl(i,j)<=vofeps) then
              phi(i,j)=dble(band_width+1)*dsqrt(PGrid%dx(i,j)**2.d0+           &
                                                PGrid%dy(i,j)**2.d0)
            elseif(vfl(i,j)>=(1.d0-vofeps-vfls(i,j))) then
              phi(i,j)=-dble(band_width+1)*dsqrt(PGrid%dx(i,j)**2.d0+          &
                                                 PGrid%dy(i,j)**2.d0)
            end if
          end if
          if(vfl(i,j)<vofeps.and.phi(i,j)<-vofeps.and.vfls(i,j)<1.d0-epsi) then
            print*,i,j
            print*,vfl(i,j),fix(i,j)
            print*,phi(i,j),phiaux(i,j)
            pause 'clsvof 1095'
          end if
        end do
      end do
    end subroutine

    subroutine Normal_Vector_Irre(PGrid,i,j,nx,ny)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid
      INTEGER,INTENT(IN):: i,j
      real(dp),INTENT(OUT):: nx,ny
      real(dp),PARAMETER:: eta = 0.075d0
      INTEGER:: dx,dy,Imv,ii
      real(dp):: qi(-1:1),qj(-1:1),nuy
      real(dp):: vx,vy,vx1,vx2,vy1,vy2,qv1,qv2,V1(2),V2(2),MaxVect,CosPar
      nuy = eta/8.d0
       ! define qi for Dx
      if(vfls(i,j)>=epsi) then
        if(vfls(i+1,j)<1.d0-epsi.and.vfls(i-1,j)<1.d0-epsi) then
          nx=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
        elseif(vfls(i+1,j)>=1.d0-epsi) then
          nx=(phi(i,j)-phi(i-1,j))/(PGrid%x(i,j)-PGrid%x(i-1,j))
        else
          nx=(phi(i+1,j)-phi(i,j))/(PGrid%x(i+1,j)-PGrid%x(i,j))
        end if
        if(vfls(i,j+1)<1.d0-epsi.and.vfls(i,j-1)<1.d0-epsi) then
          ny=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
        elseif(vfls(i,j+1)>1.d0-epsi) then
          ny=(phi(i,j)-phi(i,j-1))/(PGrid%y(i,j)-PGrid%y(i,j-1))
        else
          ny=(phi(i,j+1)-phi(i,j))/(PGrid%y(i,j+1)-PGrid%y(i,j))
        end if
        return
      end if
      if(i>=3) then
       ! for q(-1)
        vx=(phi(i,j)-phi(i-2,j))/(PGrid%x(i,j)-PGrid%x(i-2,j))
      else
        vx=(phi(i,j)-phi(i-1,j))/(PGrid%x(i,j)-PGrid%x(i-1,j))
      end if
      vy=(phi(i-1,j+1)-phi(i-1,j-1))/(PGrid%y(i-1,j+1)-PGrid%y(i-1,j-1))
      qi(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0))
       ! for q(0)
      vx=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
      vy=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
      qi(0) = dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0))
       ! for q(1)
      if(i<=ISize-2) then
        vx=(phi(i+2,j)-phi(i,j))/(PGrid%x(i+2,j)-PGrid%x(i,j))
      else
        vx=(phi(i+1,j)-phi(i,j))/(PGrid%x(i+1,j)-PGrid%x(i,j))
      end if
      vy=(phi(i+1,j+1)-phi(i+1,j-1))/(PGrid%y(i+1,j+1)-PGrid%y(i+1,j-1))
      qi(1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0))
      if(qi(-1)<eta.and.qi(1)>=eta) then
        dx=-1
  !     pause 'normal_vector_irre 790'
      elseif(qi(-1)>=eta.and.qi(1)<eta) then
        dx=1
      elseif(qi(-1)<eta.and.qi(0)<eta.and.qi(1)<eta) then
        dx=0
      elseif(qi(-1)>=eta.and.qi(0)>=eta.and.qi(1)>=eta) then
        dx=0
      else
        vx=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
        vy=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
        MaxVect=0.d0
        Imv=0
   !    Find the most parallel vector to  normal vector
        do ii = 1,8
          CosPar=dabs(Vset(ii,1)*vx+Vset(ii,2)*vy)/(dsqrt(Vset(ii,1)**2.d0+  &
                      Vset(ii,2)**2.d0)*dsqrt(vx**2.d0+vy**2.d0))
          if(CosPar>=MaxVect) then
            MaxVect=CosPar
            Imv = ii
          end if
        end do
        if(Imv==0) then
          print*,i,j
          print*,(PGrid%x(i+1,j)-PGrid%x(i-1,j))
          print*,(PGrid%y(i,j+1)-PGrid%y(i,j-1))
          print*,(phi(i+1,j)-phi(i-1,j))
          print*,(phi(i,j+1)-phi(i,j-1))
          print*,'NormalCalculation_ClsVof_1860'
          Do ii = 1,8
            CosPar=dabs(Vset(ii,1)*vx+Vset(ii,2)*vy)/(dsqrt(Vset(ii,1)**2.d0+  &
                              Vset(ii,2)**2.d0)*dsqrt(vx**2.d0+vy**2.d0))
            print*,vx,vy
            print*,
            if(CosPar>=MaxVect) then
              MaxVect=CosPar
              Imv = ii
            end if
          end do
        end if
        V1(1)=-Vset(Imv,2)
        V1(2)=Vset(Imv,1)
        V2(1)=-V1(1)
        V2(2)=-V1(2)
        if(V1(1)==-1) then
          vx1=(phi(i,j)-phi(i-1,j))/(PGrid%x(i,j)-PGrid%x(i-1,j))
          vx2=(phi(i+1,j)-phi(i,j))/(PGrid%x(i+1,j)-PGrid%x(i,j))
        elseIf(V1(1)==0) then
          vx1=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
          vx2=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
        else
          vx1=(phi(i+1,j)-phi(i,j))/(PGrid%x(i+1,j)-PGrid%x(i,j))
          vx2=(phi(i,j)-phi(i-1,j))/(PGrid%x(i,j)-PGrid%x(i-1,j))
        end if
        if(V1(2)==-1) then
          vy1=(phi(i,j)-phi(i,j-1))/(PGrid%y(i,j)-PGrid%y(i,j-1))
          vy2=(phi(i,j+1)-phi(i,j))/(PGrid%y(i,j+1)-PGrid%y(i,j))
        elseIf(V1(2)==0) then
          vy1=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
          vy2=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
        else
          vy1=(phi(i,j+1)-phi(i,j))/(PGrid%y(i,j+1)-PGrid%y(i,j))
          vy2=(phi(i,j)-phi(i,j-1))/(PGrid%y(i,j)-PGrid%y(i,j-1))
        end if
        qv1=dabs(1.d0-dsqrt(vx1**2.d0+vy1**2.d0))
        qv2=dabs(1.d0-dsqrt(vx2**2.d0+vy2**2.d0))
        if(qv1<qv2+nuy) then
          nx=vx1
          ny=vy1
          return
        else
          nx=vx2
          ny=vy2
          return
        end if
      end if
       ! define qj for Dy
      if(j>=3) then
        vy=(phi(i,j)-phi(i,j-2))/(PGrid%y(i,j)-PGrid%y(i,j-2))
      else
        vy=(phi(i,j)-phi(i,j-1))/(PGrid%y(i,j)-PGrid%y(i,j-1))
      end if
      vx=(phi(i+1,j-1)-phi(i-1,j-1))/(PGrid%x(i+1,j-1)-PGrid%x(i-1,j-1))
      qj(-1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0))
       ! define qj(0)
      vx=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
      vy=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
      qj(0)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0))
       ! define qj(1)
      vx=(phi(i+1,j+1)-phi(i-1,j+1))/(PGrid%x(i+1,j+1)-PGrid%x(i-1,j+1))
      if(j<=JSize-2) then
        vy=(phi(i,j+2)-phi(i,j))/(PGrid%y(i,j+2)-PGrid%y(i,j))
      else
        vy=(phi(i,j+1)-phi(i,j))/(PGrid%y(i,j+1)-PGrid%y(i,j))
      end if
      qj(1)=dabs(1.d0-dsqrt(vx**2.d0+vy**2.d0))
      if(qj(-1)<eta.and.qj(1)>=eta) then
        dy=-1
      elseif(qj(-1)>=eta.and.qj(1)<eta) then
        dy=1
      elseif(qj(-1)<eta.and.qj(0)<eta.and.qj(1)<eta) then
        dy=0
      elseif(qj(-1)>=eta.and.qj(0)>=eta.and.qj(1)>=eta) then
        dy=0
      else
        vx=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
        vy=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
        MaxVect=0.d0
        imv=1
      ! Find the most parallel vector to  normal vector
        do ii = 1,8
          CosPar=dabs(Vset(ii,1)*vx+Vset(ii,2)*vy)/(dsqrt(Vset(ii,1)**2.d0+  &
                              Vset(ii,2)**2.d0)*dsqrt(vx**2.d0+vy**2.d0))
          if(CosPar>=MaxVect) then
            MaxVect=CosPar
            Imv = ii
          end if
        end do
        V1(1)=-Vset(Imv,2)
        V1(2)=Vset(Imv,1)
        V2(1)=-V1(1)
        V2(2)=-V1(2)
        if(V1(1)==-1) then
          vx1=(phi(i,j)-phi(i-1,j))/(PGrid%x(i,j)-PGrid%x(i-1,j))
          vx2=(phi(i+1,j)-phi(i,j))/(PGrid%x(i+1,j)-PGrid%x(i,j))
        elseif(V1(1)==0) then
          vx1=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
          vx2=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
        else
          vx1=(phi(i+1,j)-phi(i,j))/(PGrid%x(i+1,j)-PGrid%x(i,j))
          vx2=(phi(i,j)-phi(i-1,j))/(PGrid%x(i,j)-PGrid%x(i-1,j))
        end if
        if(V1(2)==-1) then
          vy1=(phi(i,j)-phi(i,j-1))/(PGrid%y(i,j)-PGrid%y(i,j-1))
          vy2=(phi(i,j+1)-phi(i,j))/(PGrid%y(i,j+1)-PGrid%y(i,j))
        elseif(V1(2)==0) then
          vy1=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
          vy2=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
        else
          vy1=(phi(i,j+1)-phi(i,j))/(PGrid%y(i,j+1)-PGrid%y(i,j))
          vy2=(phi(i,j)-phi(i,j-1))/(PGrid%y(i,j)-PGrid%y(i,j-1))
        end if
        qv1=dabs(1.d0-dsqrt(vx1**2.d0+vy1**2.d0))
        qv2=dabs(1.d0-dsqrt(vx2**2.d0+vy2**2.d0))
        if(qv1<qv2+nuy) then
          nx=vx1
          ny=vy1
          return
        else
          nx=vx2
          ny=vy2
          return
        end if
      end if
      if(dx==-1) nx=(phi(i,j)-phi(i-1,j))/(PGrid%x(i,j)-PGrid%x(i-1,j))
      if(dx==1) nx=(phi(i+1,j)-phi(i,j))/(PGrid%x(i+1,j)-PGrid%x(i,j))
      if(dx==0) nx=(phi(i+1,j)-phi(i-1,j))/(PGrid%x(i+1,j)-PGrid%x(i-1,j))
      if(dy==-1) ny=(phi(i,j)-phi(i,j-1))/(PGrid%y(i,j)-PGrid%y(i,j-1))
      if(dy==1) ny=(phi(i,j+1)-phi(i,j))/(PGrid%y(i,j+1)-PGrid%y(i,j))
      if(dy==0) ny=(phi(i,j+1)-phi(i,j-1))/(PGrid%y(i,j+1)-PGrid%y(i,j-1))
    end subroutine Normal_Vector_Irre

    function isinsiderect(xp,yp,del_x,del_y) result(logic)
       IMPLICIT NONE
       REAL(KIND=dp),INTENT(IN):: xp,yp,del_x,del_y
       real(dp):: x1,y1,x2,y2,x3,y3,x4,y4
       real(dp):: b1,b2,b3,b4,u1,u2,u3,u4,area,areasum
       real(dp):: area1,area2,area3,area4
       logical:: logic
       x1 = -0.5*del_x
       y1 = 0.5*del_y
       x2 = 0.5*del_x
       y2 = 0.5*del_y
       x3 = 0.5*del_x
       y3 = -0.5*del_y
       x4 = -0.5*del_x
       y4 = -0.5*del_y
       b1 = dsqrt((x1-xp)**2.d0+(y1-yp)**2.d0)
       b2 = dsqrt((x2-xp)**2.d0+(y2-yp)**2.d0)
       b3 = dsqrt((x3-xp)**2.d0+(y3-yp)**2.d0)
       b4 = dsqrt((x4-xp)**2.d0+(y4-yp)**2.d0)
       u1 = 0.5d0*(del_x+b1+b2)
       u2 = 0.5d0*(del_y+b2+b3)
       u3 = 0.5d0*(del_x+b3+b4)
       u4 = 0.5d0*(del_y+b4+b1)
       area = del_x*del_y
       area1 = dsqrt(u1*(u1-del_x)*(u1-b1)*(u1-b2))
       area2 = dsqrt(u2*(u2-del_y)*(u2-b2)*(u2-b3))
       area3 = dsqrt(u3*(u3-del_x)*(u3-b3)*(u3-b4))
       area4 = dsqrt(u4*(u4-del_y)*(u4-b4)*(u4-b1))
       areasum = area1+area2+area3+area4
       if(dabs(areasum-area)<1.d-10) then
          logic = .true.
       else
          logic = .false.
       end if
    end function

    subroutine DistanceFluidCalculate(nxss,nyss,phiss,nxf,nyf,volf,dx,dy,phif)
        IMPLICIT NONE
        REAL(KIND=dp),INTENT(IN):: nxss,nyss,phiss,nxf,nyf,dx,dy,volf
        REAL(KIND=dp),INTENT(INOUT):: phif
        REAL(KIND=dp),DIMENSION(:,:):: node(6,2)
        INTEGER*4:: temp,k
        REAL(KIND=dp):: dx2,dy2,dpt(4),epsil,dt23
        REAL(KIND=dp):: MinVer(2),MaxVer(2),Mindis,Maxdis
        epsil=1.d-10
        dx2=dx/2.d0
        dy2=dy/2.d0
        temp=1
        node=0.d0
        dpt(1)=dx2*nxss-dy2*nyss+phiss
        dpt(2)=-dx2*nxss-dy2*nyss+phiss
        dpt(3)=-dx2*nxss+dy2*nyss+phiss
        dpt(4)=dx2*nxss+dy2*nyss+phiss
        if(dpt(1)>=0.d0) then
          node(temp,1)=dx2
          node(temp,2)=-dy2
          temp=temp+1
        end if
        if(dpt(1)*dpt(2)<0.d0) then
          node(temp,1)=(dy2*nyss-phiss)/nxss
          node(temp,2)=-dy2
          temp=temp+1
        end if
        if(dpt(2)>=0.d0) then
          node(temp,1)=-dx2
          node(temp,2)=-dy2
          temp=temp+1
        end if
        if(dpt(2)*dpt(3)<0.d0) then
          node(temp,1)=-dx2
          node(temp,2)=(dx2*nxss-phiss)/nyss
          temp=temp+1
        end if
        if(dpt(3)>=0.d0) then
          node(temp,1)=-dx2
          node(temp,2)=dy2
          temp=temp+1
        end if
        if(dpt(3)*dpt(4)<0.d0) then
          node(temp,1)=(-dy2*nyss-phiss)/nxss
          node(temp,2)=dy2
          temp=temp+1
        end if
        if(dpt(4)>=0.d0) then
          node(temp,1)=dx2
          node(temp,2)=dy2
          temp=temp+1
        end if
        if(dpt(4)*dpt(1)<0.d0) then
          node(temp,1)=dx2
          node(temp,2)=(-dx2*nxss-phiss)/nyss
          temp=temp+1
        end if
        node(temp,1)=node(1,1)
        node(temp,2)=node(1,2)
        Maxdis=-2.d0*dsqrt(dx**2.d0+dy**2.d0)
        Mindis=2.d0*dsqrt(dx**2.d0+dy**2.d0)
        do k=1,temp-1
          dt23=nxf*node(k,1)+nyf*node(k,2)
          if(dt23>Maxdis)then
            Maxdis=dt23
            MaxVer(1)=node(k,1)
            MaxVer(2)=node(k,2)
          end if
          if(dt23<Mindis)then
            Mindis=dt23
            MinVer(1)=node(k,1)
            MinVer(2)=node(k,2)
          end if
        end do
        MaxVer(1)=MaxVer(1)+MaxVer(1)*nxf*epsi
        MaxVer(2)=MaxVer(2)+MaxVer(2)*nyf*epsi
        MinVer(1)=MinVer(1)-MinVer(1)*nxf*epsi
        MinVer(2)=MinVer(2)-MinVer(2)*nyf*epsi
        if(temp==1.or.dabs(Node(1,1))<1.d-14.and.dabs(Node(1,2))<1.d-14) then
          print*,'problem 1266'
          print*,temp
          print*,nxss,nyss,phiss
          pause 'fuck you bugs'
        end if
        call zriddr(nxf,nyf,volf,dx,dy,phif,node,MinVer,MaxVer,1.d-14,1000)
    end subroutine DistanceFluidCalculate

    subroutine zriddr(nxf,nyf,volf,dx,dy,phif,node,x1,x2,tol,ITMAX)
       REAL(KIND=dp),INTENT(IN):: nxf,nyf,volf,dx,dy,node(6,2),x1(2),x2(2),tol
       REAL(KIND=dp),INTENT(INOUT):: phif
       INTEGER(kind=it4b),INTENT(IN):: ITMAX
       INTEGER*4:: iter
       REAL(KIND=dp):: ans(2),fh,fl,fm,fnew,s,xh(2),xl(2),xm(2),xnew(2),UNUSED
       UNUSED=-1.11d-30
       call LiqFuncVol(nxf,nyf,volf,dx,dy,node,x1,fl)
       call LiqFuncVol(nxf,nyf,volf,dx,dy,node,x2,fh)
       if((fl>tol.and.fh<-tol).or.(fl<-tol.and.fh>tol))then
         xl(:)=x1(:)
         xh(:)=x2(:)
         ans(:)=UNUSED
         do iter=1,ITMAX
           xm(:)=0.5d0*(xl(:)+xh(:))
       ! first of two function evaluation per iteration
           call LiqFuncVol(nxf,nyf,volf,dx,dy,node,xm,fm)
           s=dsqrt(fm*fm-fl*fh)
           if(s==0.d0) pause 'fuck you zriddr method'
           xnew(:)=xm(:)+(xm(:)-xl(:))*dsign(1.d0,fl-fh)*fm/s
           if(dsqrt((xnew(1)-ans(1))**2.d0+(xnew(2)-ans(2))**2.d0)<=tol) then
             phif=-(nxf*ans(1)+nyf*ans(2))
             return
           end if
           ans(:)=xnew(:)
           call LiqFuncVol(nxf,nyf,volf,dx,dy,node,xnew,fnew)
           if(dabs(fnew)<=tol) then
             phif=-(nxf*ans(1)+nyf*ans(2))
             return
           endif
           if(dsign(fm,fnew)/=fm) then
             xl(:)=xm(:)
             fl=fm
             xh(:)=ans(:)
             fh=fnew
           elseif(dsign(fl,fnew)/=fl) then
             xh(:)=ans(:)
             fh=fnew
           elseif(dsign(fh,fnew)/=fh) then
             xl=ans(:)
             fl=fnew
           else
             pause 'fuck you ridder algorithm'
           end if
           if(dsqrt((xh(1)-xl(1))**2.d0+(xh(2)-xl(2))**2.d0)<=tol) then
             phif=-(nxf*ans(1)+nyf*ans(2))
             return
           end if
         end do
         phif=-(nxf*ans(1)+nyf*ans(2))
         return
       else
         if(fl<tol.and.fl>-tol) then
           phif=-(nxf*x1(1)+nyf*x1(2))
           return
         end if
         if(fh<tol.and.fh>-tol) then
           phif=-(nxf*x2(1)+nyf*x2(2))
           return
         end if
       end if
    end subroutine zriddr

    subroutine LiqFuncVol(nxf,nyf,volf,dx,dy,node,xf,fvol)
        REAL(KIND=dp),INTENT(IN):: nxf,nyf,volf,dx,dy,node(6,2),xf(2)
        REAL(KIND=dp):: fvol,phil,nodel(7,2),pdt(6),vol
        INTEGER:: k,temp,templ
        phil=-(nxf*xf(1)+nyf*xf(2))
        temp=0
        do k = 1,6
          if(dabs(node(k,1))>1.d-14.or.dabs(node(k,2))>1.d-14) then
            temp=k
          end if
        end do
        do k = 1,temp-1
          if(k>6) then
            print*,temp
            pause 'fuck you bugs Clsvof_Mod 1468'
          end if
          pdt(k)=node(k,1)*nxf+node(k,2)*nyf+phil
        end do
        if(temp==0) then
          print*,node(1,1),node(1,2)
          print*,nxf,nyf,volf
          pause 'fuck you bugs Clsvof_Mod 1475'
        end if
        pdt(temp)=pdt(1)

    !   find the volume of liquid inside cell
        templ=1
        do k=1,temp-1
          if(pdt(k)<=0.d0) then
            nodel(templ,1)=node(k,1)
            nodel(templ,2)=node(k,2)
            templ=templ+1
          end if
          if(pdt(k)*pdt(k+1)<0.d0) then
            nodel(templ,1)=node(k,1)+(node(k+1,1)-node(k,1))*dabs(pdt(k))/     &
                                                  (dabs(pdt(k))+dabs(pdt(k+1)))
            nodel(templ,2)=node(k,2)+(node(k+1,2)-node(k,2))*dabs(pdt(k))/     &
                                                  (dabs(pdt(k))+dabs(pdt(k+1)))
            templ=templ+1
          end if
        end do
        nodel(templ,1) = nodel(1,1)
        nodel(templ,2) = nodel(1,2)
        vol = 0.d0
        do k = 1,templ-1
          vol = vol+0.5d0*(nodel(k,1)*nodel(k+1,2)-nodel(k+1,1)*nodel(k,2))
        end do
        fvol=volf-dabs(vol/(dx*dy))
    end subroutine LiqFuncVol

   ! This subroutine is used to calculate the volume fraction of liquid and solid in given cell
   ! This subroutine uses the Green theorem to calculate the volume of solid and liquid
    subroutine CellGeoCal(nxs,nys,phis,nxl,nyl,phil,dx,dy,vols,volf)
      IMPLICIT NONE
      REAL(KIND=dp),INTENT(IN):: nxs,nys,phis,nxl,nyl,phil,dx,dy
      REAL(KIND=dp),INTENT(OUT):: vols,volf
      REAL(KIND=dp),DIMENSION(:,:):: node(6,2),nodel(7,2)
      INTEGER(kind=it4b):: temp,templ,k
      REAL(KIND=dp):: dx2,dy2,pdt(6),dpt(4),vol,epsil
      epsil=1.d-24
      dx2 = dx/2.d0
      dy2 = dy/2.d0
      temp = 1
      node = 0.d0
     !  calculate the distance of all cell's nodes to interface
      dpt(1) =  dx2*nxs-dy2*nys+phis
      dpt(2) = -dx2*nxs-dy2*nys+phis
      dpt(3) = -dx2*nxs+dy2*nys+phis
      dpt(4) =  dx2*nxs+dy2*nys+phis
     !  find nodes those belong to fluid field(including the cutting nodes between
     !  interface and cell edges)
      If(dpt(1)>=0.d0) then
        node(temp,1) = dx2
        node(temp,2) = -dy2
        temp = temp+1
      End if
      If(dpt(1)*dpt(2)<0.d0) then
        node(temp,1) = (dy2*nys-phis)/nxs
        node(temp,2) = -dy2
        temp = temp+1
      End if
      If(dpt(2)>=0.d0) then
        node(temp,1) = -dx2
        node(temp,2) = -dy2
        temp = temp+1
      End if
      If(dpt(2)*dpt(3)<0.d0) then
        node(temp,1) = -dx2
        node(temp,2) = (dx2*nxs-phis)/nys
        temp = temp+1
      End if
      If(dpt(3)>=0.d0) then
        node(temp,1) = -dx2
        node(temp,2) = dy2
        temp = temp+1
      End if
      If(dpt(3)*dpt(4)<0.d0) then
        node(temp,1) = (-dy2*nys-phis)/nxs
        node(temp,2) = dy2
        temp = temp+1
      End if
      If(dpt(4)>=0.d0) then
        node(temp,1) = dx2
        node(temp,2) = dy2
        temp = temp+1
      End if
      If(dpt(4)*dpt(1)<0.d0) then
        node(temp,1) = dx2
        node(temp,2) = (-dx2*nxs-phis)/nys
        temp = temp+1
      End if
      node(temp,1) = node(1,1)
      node(temp,2) = node(1,2)
      vol = 0.d0
    ! applying Gauss theorem to find the volume of both fluid and gas
      do k = 1,temp-1
        vol = vol+0.5d0*(node(k,1)*node(k+1,2)-node(k+1,1)*node(k,2))
        pdt(k)=node(k,1)*nxl+node(k,2)*nyl+phil
      end do
    ! Volume fraction of solid in given cell
      vols = 1.d0-dabs(vol/(dx*dy))
      pdt(temp)=pdt(1)
    ! find the volume of liquid inside cell
      templ=1
      do k=1,temp-1
        if(pdt(k)<=0.d0) then
          nodel(templ,1)=node(k,1)
          nodel(templ,2)=node(k,2)
          templ=templ+1
        end if
        if(pdt(k)*pdt(k+1)<0.d0) then
          nodel(templ,1)=node(k,1)+(node(k+1,1)-node(k,1))*dabs(pdt(k))/      &
                                                 (dabs(pdt(k))+dabs(pdt(k+1)))
          nodel(templ,2)=node(k,2)+(node(k+1,2)-node(k,2))*dabs(pdt(k))/      &
                                                 (dabs(pdt(k))+dabs(pdt(k+1)))
          templ=templ+1
        end if
      end do

      nodel(templ,1) = nodel(1,1)
      nodel(templ,2) = nodel(1,2)
      vol = 0.d0
      do k = 1,templ-1
        vol = vol+0.5d0*(nodel(k,1)*nodel(k+1,2)-nodel(k+1,1)*nodel(k,2))
      end do
      volf=dabs(vol/(dx*dy))
    end subroutine CellGeoCal

    SUBROUTINE Boundary_Condition_Vof_Phi(TGrid,Time)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: TGrid
      REAL(KIND=dp),INTENT(IN):: Time
      REAL(KIND=dp),DIMENSION(:,:):: node(6,2),CutP(2,2),dpt(4)
      REAL(KIND=dp):: fx,dfx,nxx,nyy,dis,nxy,dx2,dy2,tol,vol,vos
      INTEGER(kind=it4b):: j,k,temp,templ
      do j=1,Jsize
        temp=1
        templ=1
        dx2=TGrid%dx(1,j)/2.d0
        dy2=TGrid%dy(1,j)/2.d0
        dpt(1)=TGrid%y(1,j)-dy2-Depthw-Amp0*dsin(kw*(TGrid%x(1,j)+dx2-cw0*Time))
        dpt(2)=TGrid%y(1,j)-dy2-Depthw-Amp0*dsin(kw*(TGrid%x(1,j)-dx2-cw0*Time))
        dpt(3)=TGrid%y(1,j)+dy2-Depthw-Amp0*dsin(kw*(TGrid%x(1,j)-dx2-cw0*Time))
        dpt(4)=TGrid%y(1,j)+dy2-Depthw-Amp0*dsin(kw*(TGrid%x(1,j)+dx2-cw0*Time))
        if(dpt(1)>=0.d0) then
          node(temp,1)=dx2
          node(temp,2)=-dy2
          temp=temp+1
        end if
        if(dpt(1)*dpt(2)<0.d0) then
          node(temp,1)=TGrid%x(1,j)-dx2+TGrid%dx(1,j)*dabs(dpt(2))/            &
                                                   (dabs(dpt(2))+dabs(dpt(1)))
          tol=1.d0
          do while(tol>1.d-13)
            fx=TGrid%y(1,j)-dy2-(Depthw+Amp0*dsin(kw*(node(temp,1)-cw0*Time)))
            dfx=-Amp0*kw*dcos(kw*(node(temp,1)-cw0*Time))
            tol=dabs(fx/dfx)
            node(temp,1)=node(temp,1)-fx/dfx
          end do
          node(temp,1)=node(temp,1)-TGrid%x(1,j)
          node(temp,2)=-dy2
          CutP(templ,1)=node(temp,1)
          CutP(templ,2)=node(temp,2)
          templ=templ+1
          temp=temp+1
        end if
        if(dpt(2)>=0.d0) then
          node(temp,1)=-dx2
          node(temp,2)=-dy2
          temp=temp+1
        end if
        if(dpt(2)*dpt(3)<0.d0) then
          node(temp,1)=-dx2
          node(temp,2)=Depthw+Amp0*dsin(kw*(node(temp,1)+TGrid%x(1,j)-cw0*Time))-&
                                                         TGrid%y(1,j)
          CutP(templ,1)=node(temp,1)
          CutP(templ,2)=node(temp,2)
          templ=templ+1
          temp=temp+1
        end if
        if(dpt(3)>=0.d0) then
          node(temp,1)=-dx2
          node(temp,2)=dy2
          temp=temp+1
        end if
        if(dpt(3)*dpt(4)<0.d0) then
          node(temp,1)=TGrid%x(1,j)-dx2+TGrid%dx(1,j)*dabs(dpt(3))/          &
                                       (dabs(dpt(3))+dabs(dpt(4)))
          tol=1.d0
          do while(tol>1.d-13)
            fx=TGrid%y(1,j)+dy2-(Depthw+Amp0*dsin(kw*(node(temp,1)-cw0*Time)))
            dfx=-Amp0*kw*dcos(kw*(node(temp,1)-cw0*Time))
            tol=dabs(fx/dfx)
            node(temp,1)=node(temp,1)-fx/dfx
          end do
          node(temp,1)=node(temp,1)-TGrid%x(1,j)
          node(temp,2)=dy2
          CutP(templ,1)=node(temp,1)
          CutP(templ,2)=node(temp,2)
          templ=templ+1
          temp=temp+1
        end if
        if(dpt(4)>=0.d0) then
          node(temp,1)=dx2
          node(temp,2)=dy2
          temp=temp+1
        end if
        if(dpt(4)*dpt(1)<0.d0) then
          node(temp,1)=dx2
          node(temp,2)=Depthw+Amp0*dsin(kw*(node(temp,1)+TGrid%x(1,j)-         &
                                                        cw0*Time))-TGrid%y(1,j)
          CutP(templ,1)=node(temp,1)
          CutP(templ,2)=node(temp,2)
          templ=templ+1
          temp=temp+1
        end if
        node(temp,1)=node(1,1)
        node(temp,2)=node(1,2)
        vol=0.d0
        do k=1,temp-1
          vol=vol+0.5d0*(node(k,1)*node(k+1,2)-node(k+1,1)*node(k,2))
        end do
        vol=1.d0-dabs(vol/(TGrid%dx(1,j)*TGrid%dy(1,j)))
        if(templ==3) then
          nxx=CutP(2,2)-CutP(1,2)
          nyy=CutP(2,1)-CutP(1,1)
          vos=-Amp0*cw0*kw*dcos(kw*(TGrid%x(1,j)+TGrid%dx(1,j)-cw0*Time))
          nxy=dsqrt(nxx**2.d0+nyy**2.d0)
          if(vos>0.d0) then
            nxx=-nxx/nxy
          else
            nxx=nxx/nxy
          end if
          nyy=dabs(nyy/nxy)
          dis=dabs(CutP(2,1)*CutP(1,2)-CutP(1,1)*CutP(2,2))/nxy*               &
                                                          dsign(1.d0,0.5d0-vol)
        else
          dis=TGrid%y(1,j)-(Depthw+Amp0*dsin(kw*(TGrid%x(1,j)-cw0*Time)))
          nxx=0.d0
          nyy=1.d0
        end if
        if(vol>vofeps.and.vol<-1.d0-vofeps) then
          print*,j
          print*,vol
          print*,dis
          print*,nxx,nyy
          print*,
        end if
        vfl(1,j)=vol
        phi(1,j)=dis
        nx(1,j)=nxx
        ny(1,j)=nyy
      end do
    END SUBROUTINE Boundary_Condition_Vof_Phi
 !******************************************************************************
 !*** Describing the movement of boom                                          *
 !******************************************************************************
    SUBROUTINE ObjectMovement(BoomCase,dtv)
      IMPLICIT NONE
      TYPE(SolidObject),INTENT(INOUT):: BoomCase
      REAL(KIND=dp),INTENT(IN):: dtv
      BoomCase%Posp%y=BoomCase%Posp%y+0.5d0*BoomCase%asy*dtv**2.d0+            &
                                            BoomCase%vs*dtv
      BoomCase%Posp%x=BoomCase%Posp%x+0.5d0*BoomCase%asx*dtv**2.d0+            &
                                            BoomCase%us*dtv
      BoomCase%us=BoomCase%us+BoomCase%asx*dtv
      BoomCase%vs=BoomCase%vs+BoomCase%asy*dtv
      BoomCase%XBar1=BoomCase%Posp%x-BoomCase%Wobj/2.d0
      BoomCase%XBar2=BoomCase%Posp%x+BoomCase%Wobj/2.d0
      BoomCase%YBar=BoomCase%Posp%y-dsqrt((BoomCase%Dobj/2.d0)**2.d0-          &
                   (BoomCase%Wobj/2.d0)**2.d0)-BoomCase%LBar
    END SUBROUTINE ObjectMovement
 !  Calculate the force acting on object
    SUBROUTINE ComputeForceObject(BoomCase,PGrid,PCell,VCell,TVar,ForceObj)
      IMPLICIT NONE
      TYPE(SolidObject),INTENT(INOUT):: BoomCase
      TYPE(Grid),INTENT(IN)          :: PGrid
      TYPE(Cell),INTENT(IN)          :: PCell,VCell
      TYPE(Variables),INTENT(IN)     :: TVar
      REAL(KIND=dp),INTENT(OUT)      :: ForceObj
      REAL(KIND=dp)                  :: Clp,Clf,eta,Area,Pr,Pres,Strf,nyfp,    &
                                        testArea,Maxforce,MaxP,MaxA,Maxny,     &
                                        Areatemp,SForce
      INTEGER(KIND=it4b)             :: i,j,ii,jj,temp
      Pr = 0.d0
      Clp = 0.d0
      Area = 0.d0
      TestArea=0.d0
      Maxforce=0.d0
    ! Calculate the Cdp
    ! For upper half of cylinder
      do i = ibeg,ibeg+Isize-1
        do j = jbeg,jbeg+Jsize-1
          if(PCell%vofS(i,j)<1.d0-epsi.and.PCell%vofS(i,j)>epsi) then
            SForce=0.d0
            temp=0
            Areatemp=0.d0
            if(PCell%vofS(i-1,j)>epsi.and.PCell%vofS(i-1,j)<1.d0-epsi) then
              eta=dabs(PCell%phiS(i-1,j)/(PCell%phiS(i,j)+TolDeno))
              Pres=(TVar%p(i-1,j)+TVar%p(i,j)*eta)/(eta+1.d0)
              Area=PCell%WlLh(i,j)
              nyfp=PCell%nyS(i,j)!0.5d0*(PCell%nyS(i,j)+PCell%nyS(i-1,j))
              SForce=SForce+Pres*Area*nyfp
              Areatemp=Areatemp+Area
              temp=temp+1
            end if
            if(PCell%vofS(i+1,j)>epsi.and.PCell%vofS(i+1,j)<1.d0-epsi) then
              eta=dabs(PCell%phiS(i+1,j)/(PCell%phiS(i,j)+TolDeno))
              Pres=(TVar%p(i+1,j)+TVar%p(i,j)*eta)/(eta+1.d0)
              Area=PCell%WlLh(i,j)
              nyfp=PCell%nyS(i,j)!0.5d0*(PCell%nyS(i,j)+PCell%nyS(i+1,j))
              SForce=SForce+Pres*Area*nyfp
              Areatemp=Areatemp+Area
              temp=temp+1
            end if
            if(PCell%vofS(i,j-1)>epsi.and.PCell%vofS(i,j-1)<1.d0-epsi) then
              eta=dabs(PCell%phiS(i,j-1)/(PCell%phiS(i,j)+TolDeno))
              Pres=(TVar%p(i,j-1)+TVar%p(i,j)*eta)/(eta+1.d0)
              Area=PCell%WlLh(i,j)
              nyfp=PCell%nyS(i,j)!0.5d0*(PCell%nyS(i,j)+PCell%nyS(i,j-1))
              SForce=SForce+Pres*Area*nyfp
              Areatemp=Areatemp+Area
              temp=temp+1
            end if
            if(PCell%vofS(i,j+1)>epsi.and.PCell%vofS(i,j+1)<1.d0-epsi) then
              eta=dabs(PCell%phiS(i,j+1)/(PCell%phiS(i,j)+TolDeno))
              Pres=(TVar%p(i,j+1)+TVar%p(i,j)*eta)/(eta+1.d0)
              Area=PCell%WlLh(i,j)
              nyfp=PCell%nyS(i,j)!0.5d0*(PCell%nyS(i,j)+PCell%nyS(i,j+1))
              SForce=SForce+Pres*Area*nyfp
              Areatemp=Areatemp+Area
              temp=temp+1
            end if
      !      SForce=SForce/(dble(temp)+tolDeno)
      !      Areatemp=Areatemp/(dble(temp)+tolDeno)
            SForce=TVar%p(i,j)*PCell%WlLh(i,j)*PCell%nyS(i,j)
            Areatemp=PCell%WlLh(i,j)
      !      print*,i,j
      !      print*,SForce,Areatemp,PCell%vofS(i,j),temp
            Clp=Clp+SForce
            testArea=testArea+Areatemp
            if(temp==0) then
              print*,'problem with interpolation'
              print*,i,j
              print*,PCell%vofS(i,j)
              print*,
            end if
          end if
        end do
      end do
  !    print*,'clp:',clp,testArea
  !   Some problems with bottom cells
      do i=1,Isize
        do j=1,Jsize-1
          if(PCell%VofS(i,j)>-epsi.and.PCell%vofS(i,j)<epsi) then
            if(PCell%vofS(i,j+1)<1.d0+epsi.and.PCell%vofS(i,j+1)>1.d0-epsi) then
              Pres=TVar%p(i,j)
              Area=PGrid%dx(i,j)!PCell%WlLh(i,j)
              TestArea=TestArea+Area
              nyfp=-1.d0!PCell%nyS(i,j)
              Clp=Clp+Pres*Area*nyfp
            end if
            ! 0.6d0 is the volume of solid inside a cell containing boom bar
            if(PCell%vofS(i,j+1)>0.6d0-epsi.and.PCell%vofS(i,j+1)<0.6d0+       &
                                                        epsi.and.j<Jsize-1) then
              if(PCell%vofS(i,j+2)>0.6d0-epsi.and.PCell%vofS(i,j+2)<0.6d0+     &
                                                                      epsi) then
                print*,'This is for compute force'
                print*,'There is problem with corner of boom bar'
                Pres=TVar%p(i,j)
                Area=0.6d0*PGrid%dx(i,j)!PCell%WlLh(i,j)
                TestArea=TestArea+Area
                nyfp=-1.d0!PCell%nyS(i,j)
                Clp=Clp+Pres*Area*nyfp
              end if
            end if
          end if
        end do
      end do
  !    print*,'clp 2:',clp,TestArea
    ! Wall shear stress to calculate the clf
    ! Upper half of cylinder
      clf=0.d0
      do i = ibeg,ibeg+Isize-1
        do j = jbeg,jbeg+Jsize-1
          if(VCell%VofS(i,j)>epsi.and.VCell%VofS(i,j)<1.d0-epsi.and.       &
                                                  VCell%MoExCell(i,j)/=1) then
            Strf=(TVar%v(i,j)-BoomCase%vs)/Rey/VCell%delh(i,j)
            Area=PCell%WlLh(i,j)
            Clf=Clf+Strf*area
          end if
        end do
      end do
      ForceObj=clf-clp
    END SUBROUTINE ComputeForceObject
  !******************************************************************************
  !this subroutine used to compute solid volume fraction in nextime

    SUBROUTINE SolidVolumeFraction(TGrid,TCell,BoomCase)
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN)                   :: TGrid
      TYPE(Cell),INTENT(INOUT)                :: TCell
      TYPE(SolidObject),INTENT(IN)            :: BoomCase
      TYPE(Point)                             :: Pt(0:1,0:1)
      TYPE(Point),DIMENSION(4)                :: CutP
      INTEGER(kind=it4b)                      :: i,j,ii,jj,ctr
      REAL(KIND=dp)                           :: x1,y1,dx,dy,dis,vol
      REAL(KIND=dp)                           :: tol,dpt(4)
      REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE:: nxx,nyy,VofsOld
      REAL(KIND=dp)                           :: dx1,dy1,dx2,dy2,epsil,dr1,dr2,&
                                                 nxx1,nyy1,nxy,CylBar
      allocate(nxx(Isize,Jsize))
      allocate(nyy(Isize,Jsize))
      allocate(VofsOld(Isize,Jsize))
   !  for wave only
      epsil=1.d-24
      CylBar=dsqrt((BoomCase%Dobj/2.d0)**2.d0-(BoomCase%Wobj/2.d0)**2.d0)
      VofsOld(:,:)=TCell%vofS(:,:)
      do i=1,Isize
        do j=1,Jsize
          if(isnan(TCell%vof(i,j))) then
            print*,'test'
            print*,i,j
            print*,vofsOld(i,j),TCell%VofS(i,j)
            print*,TCell%PhiS(i,j)
            print*,nxx(i,j),nyy(i,j)
            pause 'ClsVof_2668'
          end if
     !  For boom cylinder
     !     dx=TGrid%x(i,j)-BoomCase%Posp%x
     !     dy=TGrid%y(i,j)-BoomCase%Posp%y
     !     TCell%phiS(i,j)=(dsqrt(dx**2.d0+dy**2.d0)-BoomCase%Dobj/2.d0)
     !     nxx(i,j) = dx/dsqrt(dx**2.d0+dy**2.d0)
     !     nyy(i,j) = dy/dsqrt(dx**2.d0+dy**2.d0)
          do ii = 0,1
            do jj = 0,1
              Pt(ii,jj)%x = -TGrid%dx(i,j)*(0.5d0-dble(ii))+TGrid%x(i,j)
              Pt(ii,jj)%y = -TGrid%dy(i,j)*(0.5d0-dble(jj))+TGrid%y(i,j)
            end do
          end do
          ctr = 1
          CutP(:)%x = 0.d0
          CutP(:)%y = 0.d0
          call EdgeGeoCalCyl(BoomCase,Pt(0,0),Pt(0,1),CutP,ctr)
          call EdgeGeoCalCyl(BoomCase,Pt(0,1),Pt(1,1),CutP,ctr)
          call EdgeGeoCalCyl(BoomCase,Pt(1,1),Pt(1,0),CutP,ctr)
          call EdgeGeoCalCyl(BoomCase,Pt(1,0),Pt(0,0),CutP,ctr)
          if(ctr==3) then
            dx=dabs(CutP(1)%x-CutP(2)%x)
            dy=dabs(CutP(1)%y-CutP(2)%y)
            nxx(i,j)=dsign(1.d0,TGrid%x(i,j)-BoomCase%Posp%x)*dy/             &
                                                        dsqrt(dx**2.d0+dy**2.d0)
            nyy(i,j)=dsign(1.d0,TGrid%y(i,j)-BoomCase%Posp%y)*dx/             &
                                                        dsqrt(dx**2.d0+dy**2.d0)
            dx1=CutP(1)%x-TGrid%x(i,j)
            dy1=CutP(1)%y-TGrid%y(i,j)
            dx2=CutP(2)%x-TGrid%x(i,j)
            dy2=CutP(2)%y-TGrid%y(i,j)
            dis=dabs(dx1*dy2-dx2*dy1)/dsqrt(dx**2.d0+dy**2.d0)
            x1=0.5d0*(CutP(1)%x+CutP(2)%x)
            y1=0.5d0*(CutP(1)%y+CutP(2)%y)
            dr1=dsqrt((x1-BoomCase%Posp%x)**2.d0+(y1-BoomCase%Posp%y)**2.d0)
            dr2=dsqrt((TGrid%x(i,j)-BoomCase%Posp%x)**2.d0+                    &
                      (TGrid%y(i,j)-BoomCase%Posp%y)**2.d0)
            if(dr2<dr1) dis = -dis
          else
            dx=TGrid%x(i,j)-BoomCase%Posp%x
            dy=TGrid%y(i,j)-BoomCase%Posp%y+epsil
            dis=(dsqrt(dx**2.d0+dy**2.d0)-BoomCase%Dobj/2.d0)
            nxx(i,j)=dx/dsqrt(dx**2.d0+dy**2.d0)
            nyy(i,j)=dy/dsqrt(dx**2.d0+dy**2.d0)
            dis=dsign(1.d0,dis)*dmax1(dabs(dis),dabs(nxx(i,j))*TGrid%dx(i,j)/  &
                      2.d0+dabs(nyy(i,j))*TGrid%dy(i,j)/2.d0+1.d-10)
          end if
          if(ctr==3.and.dsqrt(dx**2.d0+dy**2.d0)<1.d-12) then
            dx=TGrid%x(i,j)-BoomCase%Posp%x
            dy=TGrid%y(i,j)-BoomCase%Posp%y+epsil
            dis=(dsqrt(dx**2.d0+dy**2.d0)-(BoomCase%Dobj/2.d0))
            nxx(i,j)=dx/dsqrt(dx**2.d0+dy**2.d0)
            nyy(i,j)=dy/dsqrt(dx**2.d0+dy**2.d0)
            dis=dsign(1.d0,dis)*dmax1(dabs(dis),dabs(nxx(i,j))*TGrid%dx(i,j)   &
                /2.d0+dabs(nyy(i,j))*TGrid%dy(i,j)/2.d0+1.d-10)
          end if
          TCell%phiS(i,j)=dis
          if(i==2999.and.j==1115) then
            print*,'test cut-cell for PCell inside cylinder'
            print*,i,j
            print*,dpt(1),dpt(2)
            print*,TCell%phiS(i,j)
            print*,nxx(i,j),nyy(i,j)
            print*,dsqrt((TGrid%dx(i,j)/2.d0)**2.d0+(TGrid%dy(i,j)/2.d0)**2.d0)
            print*,
          end if
     !  For region at left side of boom
          if(TGrid%x(i,j)<=BoomCase%XBar1) then
     !  For region under boom
            if(TGrid%y(i,j)<BoomCase%YBar) then
              dx=TGrid%x(i,j)-BoomCase%XBar1
              dy=TGRid%y(i,j)-BoomCase%YBar
              TCell%PhiS(i,j)=dsqrt(dx**2.d0+dy**2.d0)
              nxx(i,j)=dx/dsqrt(dx**2.d0+dy**2.d0)
              nyy(i,j)=dy/dsqrt(dx**2.d0+dy**2.d0)
     !  For region have same altitude as Boom bar
            else if(TGrid%y(i,j)<BoomCase%YBar+BoomCase%LBar) then
              dis=BoomCase%XBar1-TGrid%x(i,j)
              if(dabs(dis)<dabs(TCell%phiS(i,j))) then
                nxx(i,j)=-1.d0
                nyy(i,j)=0.d0
                TCell%phiS(i,j)=dis
              end if
            end if
      !  For region at right side of boom
          else if(TGrid%x(i,j)>=BoomCase%XBar2) then
      !  For region under boom
            if(TGrid%y(i,j)<BoomCase%YBar) then
              dx=TGrid%x(i,j)-BoomCase%XBar2
              dy=TGRid%y(i,j)-BoomCase%YBar
              TCell%PhiS(i,j)=dsqrt(dx**2.d0+dy**2.d0)
              nxx(i,j)=dx/dsqrt(dx**2.d0+dy**2.d0)
              nyy(i,j)=dy/dsqrt(dx**2.d0+dy**2.d0)
      !  For region have same altitude as Boom bar
            else if(TGrid%y(i,j)<BoomCase%YBar+BoomCase%LBar) then
              dis=TGrid%x(i,j)-BoomCase%XBar2
              if(dabs(dis)<dabs(TCell%phiS(i,j))) then
                nxx(i,j)=1.d0
                nyy(i,j)=0.d0
                TCell%phiS(i,j)=dis
              end if
            end if
       !  For region inside boom
          else
       !  For region under boom
            if(TGrid%y(i,j)<=BoomCase%YBar) then
              TCell%PhiS(i,j)=BoomCase%YBar-TGrid%y(i,j)
              nxx(i,j)=0.d0
              nyy(i,j)=-1.d0
       !  For region inside boom bar
            elseif(TGrid%y(i,j)<BoomCase%YBar+BoomCase%LBar) then
              dpt(1)=BoomCase%YBar-TGrid%y(i,j)
              dpt(2)=BoomCase%XBar1-TGrid%x(i,j)
              dpt(3)=TGrid%x(i,j)-BoomCase%XBar2
              if(dabs(dpt(1))<=dabs(dpt(2)).and.dabs(dpt(1))<=dabs(dpt(3))) then
                TCell%PhiS(i,j)=dpt(1)
                nxx(i,j)=0.d0
                nyy(i,j)=-1.d0
              end if
              if(dabs(dpt(2))<=dabs(dpt(1)).and.dabs(dpt(2))<=dabs(dpt(3))) then
                TCell%PhiS(i,j)=dpt(2)
                nxx(i,j)=-1.d0
                nyy(i,j)=0.d0
              end if
              if(dabs(dpt(3))<=dabs(dpt(1)).and.dabs(dpt(3))<=dabs(dpt(2))) then
                TCell%PhiS(i,j)=dpt(3)
                nxx(i,j)=1.d0
                nyy(i,j)=0.d0
              end if
         !  For region inside cylinder
            elseif(TGrid%y(i,j)<BoomCase%Posp%y) then
              if(dabs(TGrid%x(i,j)-0.5d0*(BoomCase%XBar1+BoomCase%XBar2))/     &
                dabs(TGrid%y(i,j)-BoomCase%Posp%y)<BoomCase%Wobj/2.d0/         &
                dsqrt((BoomCase%Dobj/2.d0)**2.d0-(BoomCase%Wobj/2.d0)**2.d0)) then
                dpt(1)=-dsqrt((TGrid%x(i,j)-BoomCase%XBar1)**2.d0+             &
                              (TGrid%y(i,j)-BoomCase%YBar-BoomCase%LBar)**2.d0)
                dpt(2)=-dsqrt((TGrid%x(i,j)-BoomCase%XBar2)**2.d0+             &
                              (TGrid%y(i,j)-BoomCase%YBar-BoomCase%LBar)**2.d0)
                if(dabs(dpt(1))<dabs(dpt(2))) then
                  TCell%PhiS(i,j)=dpt(1)
                  nxx(i,j)=(BoomCase%XBar1-TGrid%x(i,j))/dabs(dpt(1))
                  nyy(i,j)=-(TGrid%y(i,j)-BoomCase%YBar-BoomCase%LBar)/dabs(dpt(1))
                else
                  TCell%PhiS(i,j)=dpt(2)
                  nxx(i,j)=(BoomCase%XBar2-TGrid%x(i,j))/dabs(dpt(2))
                  nyy(i,j)=-(TGrid%y(i,j)-BoomCase%YBar-BoomCase%LBar)/dabs(dpt(2))
                end if
                if(i==2999.and.j==1115) then
                  print*,'test cut-cell for PCell special case'
                  print*,i,j
                  print*,dpt(1),dpt(2)
                  print*,TCell%phiS(i,j)
                  print*,nxx(i,j),nyy(i,j)
                  print*,
                end if
              end if
            end if
          end if
          if(isnan(TCell%vof(i,j))) then
            print*,'test'
            print*,vofsOld(i,j),TCell%VofS(i,j)
            print*,TCell%PhiS(i,j)
            print*,nxx(i,j),nyy(i,j)
            pause 'ClsVof_2804'
          end if
          call frac(nxx(i,j),nyy(i,j),TCell%phiS(i,j),TGrid%dx(i,j),           &
                                                      TGrid%dy(i,j),vol)
          TCell%vofS(i,j)=vol
          if(TCell%vofS(i,j)<epsi) TCell%vofS(i,j)=0.d0
          if(TCell%vofS(i,j)>=1.d0-epsi) TCell%vofS(i,j)=1.d0
          ! Correct VofS for bottom region of boom
          if(TCell%vofS(i,j)>epsi.and.TCell%vofS(i,j)<1.d0-epsi) then
            if(TGrid%y(i,j)+TGrid%dy(i,j)/2.d0-tolp>BoomCase%YBar.and.         &
               TGrid%y(i,j)-TGrid%dy(i,j)/2.d0+tolp<BoomCase%YBar) then
              if(TGrid%x(i,j)-TGrid%dx(i,j)/2.d0+tolp<BoomCase%XBar1.and.      &
                 TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-tolp>BoomCase%XBar1) then
                CutP(1)%x=BoomCase%XBar1-TGrid%x(i,j)
                CutP(1)%y=TGrid%dy(i,j)/2.d0
                CutP(2)%x=TGrid%dx(i,j)/2.d0
                CutP(2)%y=BoomCase%YBar-TGrid%y(i,j)
                nxx1=CutP(2)%y-CutP(1)%y
                nyy1=CutP(2)%x-CutP(1)%x
                TCell%vofS(i,j)=0.5d0*dabs(nxx1*nyy1)/TGrid%dx(i,j)/TGrid%dy(i,j)
                nxy=dsqrt(nxx1**2.d0+nyy1**2.d0)
                nxx(i,j)=-dabs(nxx1/nxy)
                nyy(i,j)=-dabs(nyy1/nxy)
                TCell%phiS(i,j)=dabs(CutP(2)%x*CutP(1)%y-CutP(1)%x*CutP(2)%y)/nxy
              end if
              if(TGrid%x(i,j)-TGrid%dx(i,j)/2.d0+tolp<BoomCase%XBar2.and.      &
                 TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-tolp>BoomCase%XBar2) then
                CutP(1)%x=BoomCase%XBar2-TGrid%x(i,j)
                CutP(1)%y=TGrid%dy(i,j)/2.d0
                CutP(2)%x=-TGrid%dx(i,j)/2.d0
                CutP(2)%y=BoomCase%YBar-TGrid%y(i,j)
                nxx1=CutP(2)%y-CutP(1)%y
                nyy1=CutP(2)%x-CutP(1)%x
                TCell%vofS(i,j)=0.5d0*dabs(nxx1*nyy1)/TGrid%dx(i,j)/TGrid%dy(i,j)
                nxy=dsqrt(nxx1**2.d0+nyy1**2.d0)
                nxx(i,j)=dabs(nxx1/nxy)
                nyy(i,j)=-dabs(nyy1/nxy)
                TCell%phiS(i,j)=dabs(CutP(2)%x*CutP(1)%y-CutP(1)%x*CutP(2)%y)/nxy
              end if
            end if
          end if
          ! Correct VofS for region containing cylinder and bar
          if(TCell%vofS(i,j)>epsi.and.TCell%vofS(i,j)<1.d0-epsi) then
            if(TGrid%y(i,j)+TGrid%dy(i,j)/2.d0-tolp>BoomCase%Posp%y-CylBar.and. &
               TGrid%y(i,j)-TGrid%dy(i,j)/2.d0+tolp<BoomCase%Posp%y-CylBar)then
              if(TGrid%x(i,j)-TGrid%dx(i,j)/2.d0+tolp<BoomCase%XBar1.and.       &
                 TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-tolp>BoomCase%XBar1) then
                CutP(1)%x=BoomCase%XBar1-TGrid%x(i,j)
                CutP(1)%y=-TGrid%dy(i,j)/2.d0
                if(dsqrt((TGrid%x(i,j)-TGrid%dx(i,j)/2.d0-BoomCase%Posp%x)**2.d0+&
                  (TGrid%y(i,j)+TGrid%dy(i,j)/2.d0-BoomCase%Posp%y)**2.d0)<    &
                  BoomCase%Dobj/2.d0) then
                  CutP(2)%x=-TGrid%dx(i,j)/2.d0
                  CutP(2)%y=-dsqrt((BoomCase%Dobj/2.d0)**2.d0-                 &
                    (TGrid%x(i,j)-BoomCase%Posp%x-TGrid%dx(i,j)/2.d0)**2.d0)-  &
                    (TGrid%y(i,j)-BoomCase%Posp%y)
                  nxx1=CutP(2)%y-CutP(1)%y
                  nyy1=CutP(2)%x-CutP(1)%x
                  TCell%vofS(i,j)=1.d0-0.5d0*dabs(nxx1*nyy1)/TGrid%dx(i,j)/    &
                                                                  TGrid%dy(i,j)
                  nxy=dsqrt(nxx1**2.d0+nyy1**2.d0)
                  nxx(i,j)=-dabs(nxx1/nxy)
                  nyy(i,j)=-dabs(nyy1/nxy)
                  TCell%phiS(i,j)=dabs(CutP(2)%x*CutP(1)%y-                    &
                       CutP(1)%x*CutP(2)%y)/nxy*dsign(1.d0,0.5d0-TCell%vofS(i,j))
                else
                  CutP(2)%x=-dsqrt((BoomCase%Dobj/2.d0)**2.d0-                  &
                    (BoomCase%Posp%y-TGrid%y(i,j)-TGrid%dy(i,j)/2.d0)**2.d0)-  &
                    (TGrid%x(i,j)-BoomCase%Posp%x)
                  CutP(2)%y=TGrid%dy(i,j)/2.d0
                  nxx1=CutP(2)%y-CutP(1)%y
                  nyy1=CutP(2)%x-CutP(1)%x
                  TCell%vofS(i,j)=0.5d0*dabs(TGrid%dx(i,j)-CutP(1)%x-CutP(2)%x)&
                                                                /TGrid%dx(i,j)
                  nxy=dsqrt(nxx1**2.d0+nyy1**2.d0)
                  nxx(i,j)=-dabs(nxx1/nxy)
                  nyy(i,j)=-dabs(nyy1/nxy)
                  TCell%phiS(i,j)=dabs(CutP(2)%x*CutP(1)%y-                    &
                       CutP(1)%x*CutP(2)%y)/nxy*dsign(1.d0,0.5d0-TCell%vofS(i,j))
                end if
              end if
              if(TGrid%x(i,j)-TGrid%dx(i,j)/2.d0+tolp<BoomCase%XBar2.and.      &
                 TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-tolp>BoomCase%XBar2)then
                CutP(1)%x=BoomCase%XBar2-TGrid%x(i,j)
                CutP(1)%y=-TGrid%dy(i,j)/2.d0
                if(dsqrt((TGrid%x(i,j)+TGrid%dx(i,j)/2.d0-BoomCase%Posp%x)**2.d0+&
                  (TGrid%y(i,j)+TGrid%dy(i,j)/2.d0-BoomCase%Posp%y)**2.d0)<    &
                  BoomCase%Dobj/2.d0) then
                  CutP(2)%x=TGrid%dx(i,j)/2.d0
                  CutP(2)%y=-dsqrt((BoomCase%Dobj/2.d0)**2.d0-                 &
                    (TGrid%x(i,j)-BoomCase%Posp%x+TGrid%dx(i,j)/2.d0)**2.d0)-  &
                    (TGrid%y(i,j)-BoomCase%Posp%y)
                  nxx1=CutP(2)%y-CutP(1)%y
                  nyy1=CutP(2)%x-CutP(1)%x
                  TCell%vofS(i,j)=1.d0-0.5d0*dabs(nxx1*nyy1)/TGrid%dx(i,j)/    &
                                                             TGrid%dy(i,j)
                  nxy=dsqrt(nxx1**2.d0+nyy1**2.d0)
                  nxx(i,j)=dabs(nxx1/nxy)
                  nyy(i,j)=-dabs(nyy1/nxy)
                  TCell%phiS(i,j)=dabs(CutP(2)%x*CutP(1)%y-                    &
                       CutP(1)%x*CutP(2)%y)/nxy*dsign(1.d0,0.5d0-TCell%vofS(i,j))
                else
                  CutP(2)%x=dsqrt((BoomCase%Dobj/2.d0)**2.d0-                  &
                    (BoomCase%Posp%y-TGrid%y(i,j)-TGrid%dy(i,j)/2.d0)**2.d0)-  &
                    (TGrid%x(i,j)-BoomCase%Posp%x)
                  CutP(2)%y=TGrid%dy(i,j)/2.d0
                  nxx1=CutP(2)%y-CutP(1)%y
                  nyy1=CutP(2)%x-CutP(1)%x
                  TCell%vofS(i,j)=0.5d0*dabs(CutP(1)%x+CutP(2)%x+TGrid%dx(i,j))&
                                                                /TGrid%dx(i,j)
                  nxy=dsqrt(nxx1**2.d0+nyy1**2.d0)
                  nxx(i,j)=dabs(nxx1/nxy)
                  nyy(i,j)=-dabs(nyy1/nxy)
                  TCell%phiS(i,j)=dabs(CutP(2)%x*CutP(1)%y-                    &
                       CutP(1)%x*CutP(2)%y)/nxy*dsign(1.d0,0.5d0-TCell%vofS(i,j))
                end if
              end if
            end if
          end if
          if(i==2999.and.j==1114) then
            print*,'test inside boom'
            print*,TCell%phiS(i,j)
            print*,TCell%vofS(i,j)
            print*,TCell%nxS(i,j)
            print*,TCell%nyS(i,j)
            print*,
          end if
          TCell%nxS(i,j)=nxx(i,j)
          TCell%nyS(i,j)=nyy(i,j)
          if(i==2999.and.j==1115) then
            print*,'after inside boom'
            print*,i,j
            print*,TCell%phiS(i,j)
            print*,TCell%vofS(i,j)
            print*,TCell%nxS(i,j)
            print*,TCell%nyS(i,j)
            print*,
          end if
          if(TCell%vofS(i,j)<epsi) TCell%vofS(i,j)=0.d0
          if(TCell%vofS(i,j)>=1.d0-epsi) TCell%vofS(i,j)=1.d0
          if((TCell%Vof(i,j)+VofsOld(i,j)>=1.d0-epsi).or.                    &
                             (TCell%Vof(i,j)+TCell%VofS(i,j))>=1.d0-epsi) then
            TCell%Vof(i,j)=1.d0-TCell%VofS(i,j)
          end if
          if(VofsOld(i,j)<1.d0-epsi.and.                                       &
            TCell%Vof(i,j)+VofsOld(i,j)<1.d0-epsi.and.                       &
            TCell%Vof(i,j)+TCell%VofS(i,j)<1.d0-epsi) then
            TCell%Vof(i,j)=TCell%Vof(i,j)/(1.d0-VofsOld(i,j)+TolDeno)*         &
                                  (1.d0-TCell%VofS(i,j))
          end if
          if(isnan(TCell%vof(i,j))) then
            print*,'test'
            print*,vofsOld(i,j),TCell%VofS(i,j)
            print*,TCell%PhiS(i,j)
            print*,nxx(i,j),nyy(i,j)
            pause 'ClsVof_2875'
          end if
        end do
      end do
      deallocate(nxx)
      deallocate(nyy)
    END SUBROUTINE SolidVolumeFraction

    SUBROUTINE EdgeGeoCalCyl(BoomCase,Pt1,Pt2,CutP,ctr)
       IMPLICIT NONE
       TYPE(SolidObject),INTENT(IN)          :: BoomCase
       TYPE(point),INTENT(IN)                :: Pt1,Pt2
       TYPE(point),DIMENSION(4),INTENT(INOUT):: CutP
       INTEGER(KIND=it4b),INTENT(INOUT)      :: ctr
       REAL(KIND=dp)                         :: epsil,Lvs_Pt1,Lvs_Pt2,EdAr
       Lvs_Pt1=dsqrt((Pt1%x-BoomCase%Posp%x)**2.d0+                            &
                     (Pt1%y-BoomCase%Posp%y)**2.d0)-BoomCase%Dobj/2.d0
       Lvs_Pt2=dsqrt((Pt2%x-BoomCase%Posp%x)**2.d0+                            &
                     (Pt2%y-BoomCase%Posp%y)**2.d0)-BoomCase%Dobj/2.d0
       epsil = 1.d-20
       if(Lvs_Pt1>=epsil.and.Lvs_Pt2<epsil) then
         if(dabs(pt1%x-pt2%x)<=1.d-10) then
           EdAr=pt1%y-dsign(1.d0,pt1%y-BoomCase%Posp%y)*                       &
                dsqrt((BoomCase%Dobj/2.d0)**2.d0-(pt1%x-BoomCase%Posp%x)**2.d0)&
                -BoomCase%Posp%y
           CutP(ctr)%y=pt1%y-EdAr
           CutP(ctr)%x=pt1%x
           ctr=ctr+1
         else
           EdAr=pt1%x-dsign(1.d0,pt1%x-BoomCase%Posp%x)*                       &
                dsqrt((BoomCase%Dobj/2.d0)**2.d0-(pt1%y-BoomCase%Posp%y)**2.d0)&
                -BoomCase%Posp%x
           CutP(ctr)%x=pt1%x-EdAr
           CutP(ctr)%y=pt1%y
           ctr=ctr+1
         end if
       elseIf(Lvs_Pt1<epsil.and.Lvs_Pt2>=epsil) then
         if(dabs(pt1%x-pt2%x)<=1.d-10) then
           EdAr=pt2%y-dsign(1.d0,pt2%y-BoomCase%Posp%y)*                       &
                dsqrt((BoomCase%Dobj/2.d0)**2.d0-(pt2%x-BoomCase%Posp%x)**2.d0)&
                -BoomCase%Posp%y
           CutP(ctr)%y=pt2%y-EdAr
           CutP(ctr)%x=pt1%x
           ctr=ctr+1
         else
           EdAr=pt2%x-dsign(1.d0,pt2%x-BoomCase%Posp%x)*                       &
                dsqrt((BoomCase%Dobj/2.d0)**2.d0-(pt2%y-BoomCase%Posp%y)**2.d0)&
                -BoomCase%Posp%x
           CutP(ctr)%x=pt2%x-EdAr
           CutP(ctr)%y=pt2%y
           ctr=ctr+1
         end if
       end if
     END SUBROUTINE EdgeGeoCalCyl
End module Clsvof
