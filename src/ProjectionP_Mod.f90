Module ProjectionP
 !! Description:
 !! The module compute the projection step for pressure.
 !! Method:
 !! The linear system will be solved by HYPRE.
 ! Current Code Owner: SIMCOFlow
 ! Code Description:
 ! Language: Fortran 90.
 ! Software Standards: "European Standards for Writing and
 ! Documenting Exchangeable Fortran 90 Code".
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! Author : BjSon Tung Dang
 !        : NTNU,SINTEF
 ! Date : 20.09.2019
    USE PrecisionVar
    USE Mesh, ONLY : Grid, Cell, ibeg, jbeg, Isize, Jsize
    USE StateVariables, ONLY : TVariables
    USE Constants, ONLY : epsi, roa, row
    USE PredictorUV, ONLY : Predictor, PoissonCoefficient
    USE MPI
    use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
    IMPLICIT NONE
    PRIVATE
    REAL(KIND=dp),DIMENSION(:,:),POINTER:: p,u,v
    REAL(KIND=dp),PARAMETER:: amp = 1.d0
    REAL(KIND=dp),PARAMETER:: alp=0.4d0,bet=0.5d0
    REAL(KIND=dp),DIMENSION(:,:,:),ALLOCATABLE:: PoCoef,HJump   ! the coefficient for Poisson solving
    TYPE,PUBLIC:: Projection
      REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE:: Pp
    END TYPE
    PUBLIC:: PoissonEquationSolver
    INTERFACE PoissonEquationSolver
      MODULE PROCEDURE PoissonEquationSolver
    END INTERFACE PoissonEquationSolver
    CONTAINS
    SUBROUTINE PoissonEquationSolver(PGrid,UGrid,VGrid,PCellO,PCell,UCell,     &
                                     VCell,TVar,TPred,PU,PV,Proj,vb,dt,itt)
      !! The subroutine is solving the Poisson like equation	
      IMPLICIT NONE
      TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      TYPE(Cell),INTENT(IN):: PCellO,PCell,UCell,VCell
      TYPE(TVariables),INTENT(IN),target:: TVar
=======
      !! The grid	 	
      TYPE(Cell),INTENT(IN):: PCellO
      !! The old pressure cell	
      TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      !! The current pressure cell 	
      TYPE(Cell),INTENT(IN):: UCell,VCell
      !! The current velocity cell 
      TYPE(Variables),INTENT(IN),target:: TVar
      !! The state variables	
>>>>>>> origin/CommentedCode
      TYPE(Predictor),INTENT(IN),target:: TPred
      !! The predicted velocities
      TYPE(Projection),INTENT(INOUT),target:: Proj
      !! The projection pressure  	
      TYPE(PoissonCoefficient),INTENT(IN):: PU,PV
      !! The coefficient for Poisson equation 
      REAL(KIND=dp),INTENT(IN):: vb
      !! The boundary velocity	 
      INTEGER(kind=it8b),INTENT(IN):: itt
      !! The iteration step	
      REAL(KIND=dp),INTENT(IN):: dt
      !! The time step	
      INTEGER*8:: A,parcsr_A,b,par_b,x,par_x,solver,precond
      INTEGER(kind=it4b):: num_iterations,i,j
      REAL(KIND=dp):: final_res_norm,tol,resi
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: matr
      REAL(KIND=dp),DIMENSION(:,:),allocatable:: rhm
      allocate(matr(Isize,Jsize,5))
      allocate(rhm(Isize,Jsize))
      allocate(PoCoef(Isize,Jsize,4)) ! the order of the face: S-1;W-2,E-3,N-4
      allocate(HJump(Isize,Jsize,4))
      p => TVar%p
      u => TPred%u
      v => TPred%v
      call ComputePossionMatrixCoefficient(PGrid,UGrid,VGrid,PCell,UCell,      &
                                                             VCell,PU,PV, TVar%Roref)
      call SetBasicSolver(solver,precond)
      ! call SetBasicSolver(solver=solver,ierr=ierr) 
      ! The order of boundary condition WB,EB,SB,NB. 1 represent  
      !
      call SetPoissonMatrix(A,parcsr_A,PGrid,PCell,1,1,1,0,matr,itt)
      call SetPoissonVectors(b,x,par_b,par_x,PGrid,PCellO,PCell,vb,dt,rhm,itt)
      call HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x,ierr)
      call HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x,ierr)
      ! Run info - needed logging turned on
      call HYPRE_ParCSRPCGGetNumIterations(solver,num_iterations,ierr)
      call HYPRE_ParCSRPCGGetFinalRelative(solver,final_res_norm,ierr)
      ! Call HYPRE_ParCSRPCGGetResidual(solver,tol,ierr)
      call DeltaPressureGetValues(x,PCell,Proj)
      call DeltaPressureBoundaryCondition(Proj,1,1,1,0)
      call HYPRE_IJMatrixDestroy(A,ierr)
      call HYPRE_IJVectorDestroy(b,ierr)
      call HYPRE_IJVectorDestroy(x,ierr)
      call HYPRE_BoomerAMGDestroy(precond,ierr)
      call HYPRE_ParCSRPCGDestroy(solver,ierr)
      do i = ibeg,Isize
        do j = jbeg,Jsize
          resi = 0.d0
          if(j>1) then
            resi=resi+matr(i,j,1)*Proj%Pp(i,j-1)
          end if
          if(i>1) then
            resi=resi+matr(i,j,2)*Proj%Pp(i-1,j)
          end if
          resi=resi+matr(i,j,3)*Proj%Pp(i,j)
          if(i<Isize) then
            resi=resi+matr(i,j,4)*Proj%Pp(i+1,j)
          end if
          if(j<Jsize) then
            resi=resi+matr(i,j,5)*Proj%Pp(i,j+1)
          end if
          if(dabs(resi-rhm(i,j))>1.d-10.and.PCell%Posnu(i,j)/=-1) then
            print*,'Problem start'
            print*,resi,rhm(i,j),dabs(resi-rhm(i,j))
            print*,i,j
            print*, 'ProjectionP_86'
            read(stdin,*)
          end if
        end do
      end do
      nullify(p)
      nullify(u)
      nullify(v)
      deallocate(matr,rhm)
      deallocate(PoCoef,HJump)
    END SUBROUTINE PoissonEquationSolver

    SUBROUTINE SetBasicSolver(solver,precond)
      !! The subroutine set up the solver and preconditioning for HYPRE	
      IMPLICIT NONE
      INTEGER*8,INTENT(INOUT):: solver
      !! The Id for solver method	
      INTEGER*8,INTENT(INOUT),optional:: precond
      !! The Id for preconditioning method	
      ! Set up and use a solver
      call HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD,solver,ierr)
      ! Set some PARAMETERs
      call HYPRE_ParCSRPCGSetMaxIter(solver,50,ierr)
      ! call HYPRE_ParCSRPCGSetAbsoluteTol(solver,1.d-14,ierr)
      ! call HYPRE_ParCSRPCGSetResidualTol(solver,1.0d-14,ierr)
      call HYPRE_ParCSRPCGSetTol(solver,1.0d-20,ierr)
      call HYPRE_ParCSRPCGSetTwoNorm(solver,0,ierr)
      ! call HYPRE_ParCSRPCGSetPrintLevel(solver,2,ierr)
      call HYPRE_ParCSRPCGSetLogging(solver,1,ierr)
      ! Now set up the AMG preconditioner and specify any PARAMETERs
      if(present(precond)) then
        call HYPRE_BoomerAMGCreate(precond,ierr)
        ! Set some PARAMETERs
        ! Print less solver info since a preconditioner
        ! call HYPRE_BoomerAMGSetPrintLevel(precond,1,ierr);
        ! falgout coarsening
        call HYPRE_BoomerAMGSetCoarsenTYPE(precond,6,ierr)
        ! SYMMETRIC G-S/Jacobi hybrid relaxation
        call HYPRE_BoomerAMGSetRelaxTYPE(precond,6,ierr)
        ! Sweeeps on each level
        call HYPRE_BoomerAMGSetNumSweeps(precond,1,ierr)
        ! conv. tolerance
        call HYPRE_BoomerAMGSetTol(precond,0.0d0,ierr)
        ! do only one iteration!
        call HYPRE_BoomerAMGSetMaxIter(precond,10,ierr)
        ! set amg as the pcg preconditioner
        ! precond_id = 2
        call HYPRE_ParCSRPCGSetPrecond(solver,2,precond,ierr)
      end if
    END SUBROUTINE SetBasicSolver

    SUBROUTINE SetPoissonMatrix(A,parcsr_A,PGrid,PCell,WB,EB,SB,NB,matr,itt)
      !! The subroutine set up matrix coefficients for linear system
      IMPLICIT NONE
      INTEGER*8,INTENT(INOUT):: A,parcsr_A
      !! The Id for matrix
      TYPE(Grid),INTENT(IN):: PGrid
      !! The grid	
      TYPE(Cell),INTENT(IN):: PCell
      !! The cell
      INTEGER(kind=it4b),INTENT(IN):: WB,EB,SB,NB 
      !! Boundary condition for west face, east face,
      !! south face, north face, 0: Dirichlet, 1: Neumann
      INTEGER(kind=it8b),INTENT(IN):: itt
      !! The iteration steps
      INTEGER(kind=it4b):: nnz,ictr,ilower,iupper,cols(0:4)
      INTEGER(kind=it4b):: i,j,ii,jj
      REAL(KIND=dp):: values(0:4)
      REAL(KIND=dp):: dx,dy,test,nesu,diag,tol
      REAL(KIND=dp),DIMENSION(:,:,:),allocatable:: matr
      ilower = 0
      iupper = PCell%ExtCell
      tol = 1.d-24
      ! Create and Set up matrix
      call HYPRE_IJMatrixCreate(MPI_COMM_WORLD,ilower,iupper,ilower,iupper,    &
                                                                        A,ierr)
      call HYPRE_IJMatrixSetObjectTYPE(A,HYPRE_PARCSR,ierr)
      call HYPRE_IJMatrixInitialize(A,ierr)
      ! Now go through my local rows and set the matrix entries.
      ! Each row has at most 5 entries. For example, if n=3:
      !
      ! A = [M -I 0; -I M -I; 0 -I M]
      ! M = [4 -1 0; -1 4 -1; 0 -1 4]
      !
      ! Note that here we are setting one row at a time, though
      ! one could set all the rows together (see the User's Manual).
      ! mindiag=1.d10
      do i = ibeg,ibeg+Isize-1
        do j = jbeg,jbeg+Jsize-1
          if(PCell%Posnu(i,j)/=-1) then
            nesu=0.d0
            diag=0.d0
            matr(i,j,:)=0.d0
            dx=PGrid%dx(i,j)
            dy=PGrid%dy(i,j)
            ictr=Pcell%Posnu(i,j)
            nnz=0
            values=0.d0
            cols=0
            ! Bottom of current cell
            if(j>jbeg) then
              if(PCell%Posnu(i,j-1)/=-1) then
                cols(nnz)=PCell%Posnu(i,j-1)
                values(nnz)=-PoCoef(i,j,1)*PCell%SEdge_Area(i,j)
                matr(i,j,1)=values(nnz)
                nesu=nesu+dabs(values(nnz))
                nnz=nnz+1
              end if
            end if
            ! West of current cell
            if(i>ibeg) then
              if(PCell%Posnu(i-1,j)/=-1) then
                cols(nnz)=PCell%Posnu(i-1,j)
                values(nnz)=-PoCoef(i,j,2)*PCell%WEdge_Area(i,j)
                nesu=nesu+dabs(values(nnz))
                matr(i,j,2)=values(nnz)
                nnz=nnz+1
              end if
            end if
            ! Set the diagonal cell
            cols(nnz)=PCell%Posnu(i,j)
            values(nnz)=PoCoef(i,j,1)*PCell%SEdge_Area(i,j)+                   &
                        PoCoef(i,j,2)*PCell%WEdge_Area(i,j)+                   &
                        PoCoef(i,j,3)*PCell%EEdge_Area(i,j)+                   &
                        PoCoef(i,j,4)*PCell%NEdge_Area(i,j)
            if(isnan(Values(nnz)).or.dabs(Values(nnz))>1.d10) then
              print*,i,j,Values(nnz)
              print*,PoCoef(i,j,1),PoCoef(i,j,2),PoCoef(i,j,3),PoCoef(i,j,4)
              print*, 'ProjectP_195 Bugs, you are again!'
              read(stdin,*)
            end if
            ! Apply boundary condition for matrix
            if(SB==1.and.j==jbeg) then
              values(nnz)=values(nnz)-PoCoef(i,j,1)*PCell%SEdge_Area(i,j)
            end if
            if(WB==1.and.i==ibeg) then
              values(nnz)=values(nnz)-PoCoef(i,j,2)*PCell%WEdge_Area(i,j)
            end if
            if(EB==1.and.i==ibeg+Isize-1) then
              values(nnz)=values(nnz)-PoCoef(i,j,3)*PCell%EEdge_Area(i,j)
            end if
            if(NB==1.and.j==jbeg+Jsize-1) then
              values(nnz)=values(nnz)-PoCoef(i,j,4)*PCell%NEdge_Area(i,j)
            end if
            values(nnz)=values(nnz)+dsign(1.d0,values(nnz))*tol
            diag=dabs(values(nnz))
            matr(i,j,3)=values(nnz)
            nnz=nnz+1
            ! East of current cell
            if(i<ibeg+Isize-1) then
              if(PCell%Posnu(i+1,j)/=-1) then
                cols(nnz)=PCell%Posnu(i+1,j)
                values(nnz)=-PoCoef(i,j,3)*PCell%EEdge_Area(i,j)
                nesu=nesu+dabs(values(nnz))
                matr(i,j,4)=values(nnz)
                nnz=nnz+1
              end if
            end if
            ! North of current cell
            if(j<jbeg+Jsize-1) then
              if(PCell%Posnu(i,j+1)/=-1) then
                cols(nnz)=PCell%Posnu(i,j+1)
                values(nnz)=-PoCoef(i,j,4)*PCell%NEdge_Area(i,j)
                nesu=nesu+dabs(values(nnz))
                matr(i,j,5)=values(nnz)
                nnz=nnz+1
              end if
            end if
            call HYPRE_IJMatrixSetValues(A,1,nnz,ictr,cols,values,ierr)
          end if
        end do
      end do
      call HYPRE_IJMatrixAssemble(A,ierr)
      call HYPRE_IJMatrixGetObject(A,parcsr_A,ierr)
    end subroutine SetPoissonMatrix

    subroutine SetPoissonVectors(b,x,par_b,par_x,PGrid,PCellO,PCell,vb,dt,rhm,itt)
      !! The subroutine setup the right hand side vector and roots for linear system
        INTEGER*8,INTENT(IN):: b,par_b
        !! The id of right hand side vector
        INTEGER*8,INTENT(IN):: x,par_x
	!! The id of roots
        TYPE(Grid),INTENT(IN):: PGrid
	!! The grid
        TYPE(Cell),INTENT(IN):: PCellO
	!! The old pressure cell 
	TYPE(Cell),INTENT(IN):: PCell
        !! The current pressure cell         
	REAL(KIND=dp),INTENT(IN):: vb
	!! The boundary velocity 
        REAL(KIND=dp),INTENT(IN):: dt
	!! The time step size
        INTEGER(kind=it8b),INTENT(IN):: itt
        !! The iteration steps
        INTEGER(kind=it4b):: i,j,ii,jj
        INTEGER:: ilower,iupper,ictr,local_size
        REAL(KIND=dp):: dx,dy,beta(2),maxvect
        INTEGER(kind=it4b),DIMENSION(:),allocatable:: rows
        REAL(KIND=dp),DIMENSION(:),allocatable:: rhs,xval
        REAL(KIND=dp),DIMENSION(:,:),allocatable:: ExtFlux,rhm
        ilower = 0
        iupper = PCell%ExtCell
        local_size = iupper-ilower+1 ! the number of rows
        ! In here, we apply boundary condition for deltaP with its values is 0 at
        ! all boundary. therefore, we do not need to set boundary in vector b
        allocate(rhs(0:PCell%ExtCell))
        allocate(xval(0:PCell%ExtCell))
        allocate(rows(0:PCell%ExtCell))
        allocate(ExtFlux(ibeg:ibeg+Isize-1,jbeg+jbeg+Jsize-1))
        call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,b,ierr)
        call HYPRE_IJVectorSetObjectTYPE(b,HYPRE_PARCSR,ierr)
        call HYPRE_IJVectorInitialize(b,ierr)

        call HYPRE_IJVectorCreate(MPI_COMM_WORLD,ilower,iupper,x,ierr)
        call HYPRE_IJVectorSetObjectTYPE(x,HYPRE_PARCSR,ierr)
        call HYPRE_IJVectorInitialize(x,ierr)
        ExtFlux(:,:) = 0.d0
        ! maxvect=0.d0
        do i = ibeg,ibeg+Isize-1
          do j = jbeg,jbeg+Jsize-1
            ictr=PCell%PosNu(i,j)
            if(ictr/=-1) then
              rhm(i,j) = 0.d0
              dx=Pgrid%dx(i,j)
              dy=Pgrid%dy(i,j)
              rhs(ictr)=-(u(i,j)*PCell%EEdge_Area(i,j)-                        &
                          u(i-1,j)*PCell%WEdge_Area(i,j))*dy-                  &
                         (v(i,j)*PCell%NEdge_Area(i,j)-                        &
                          v(i,j-1)*PCell%SEdge_Area(i,j))*dx                   &
                         +vb*PCell%nyS(i,j)*PCell%WlLh(i,j)
                      ! -((1.d0-PCell%vofS(i,j))-(1.d0-PCellO%vofS(i,j)))*dx*dy
              xval(ictr)=0.d0
              rows(ictr)=ilower+ictr
              rhm(i,j)=rhs(ictr)
              if(isnan(rhs(ictr)).or.dabs(rhs(ictr))>1.d5) then
                print*,u(i,j),u(i-1,j),v(i,j),v(i,j-1),i,j
                print*, 'bugs, projection 476'
                read(stdin,*)
              endif
            endif
          end do
        end do
        call HYPRE_IJVectorSetValues(b,local_size,rows,rhs,ierr)
        call HYPRE_IJVectorSetValues(x,local_size,rows,xval,ierr)
        call HYPRE_IJVectorAssemble(b,ierr)
        call HYPRE_IJVectorAssemble(x,ierr)
        ! get the x and b objects
        call HYPRE_IJVectorGetObject(b,par_b,ierr)
        call HYPRE_IJVectorGetObject(x,par_x,ierr)
        deallocate(rhs)
        deallocate(xval)
        deallocate(rows)
        deallocate(ExtFlux)
    end subroutine SetPoissonVectors

    subroutine DeltaPressureGetValues(x,PCell,Projp)
      !! The subroutine will get the roots from HYPRE. 	
        INTEGER*8,INTENT(IN):: x
        TYPE(Cell),INTENT(IN):: PCell
        TYPE(Projection),INTENT(INOUT):: Projp
        INTEGER(kind=it4b):: i,j
        INTEGER(kind=it4b):: ilower,iupper,local_size,ctr
        INTEGER(kind=it4b),DIMENSION(:),allocatable:: rows
        REAL(KIND=dp),DIMENSION(:),allocatable:: values
        ilower = 0
        iupper = PCell%ExtCell
        local_size = PCell%ExtCell+1 ! number of element
        allocate(values(ilower:iupper))
        allocate(rows(ilower:iupper))
        do i = ibeg,ibeg+Isize-1
          do j = jbeg,jbeg+Jsize-1
            if(PCell%PosNu(i,j)/=-1) then
              rows(PCell%PosNu(i,j)) = PCell%PosNu(i,j)+ilower
            end if
          end do
        end do
        call HYPRE_IJVectorGetValues(x,local_size,rows,values,ierr)
        ctr = 0
        do i = ibeg,ibeg+Isize-1
          do j = jbeg,jbeg+Jsize-1
            if(PCell%PosNu(i,j)==ctr) then
              Projp%Pp(i,j) = values(ctr)
              ctr = ctr+1
              if(isnan(Projp%Pp(i,j)).or.projp%Pp(i,j)>1.d10) then
                print*,Projp%Pp(i,j)
                print*, 'ProjectionP_Mod 356'
                read(stdin,*)
              end if
            else
              Projp%Pp(i,j) = 0.d0
            end if
          end do
        end do
        deallocate(values,rows)
    end subroutine DeltaPressureGetValues

    subroutine ComputePossionMatrixCoefficient(PGrid,UGrid,VGrid,PCell,UCell,  &
<<<<<<< HEAD
                                                                 VCell,PU,PV, Roref)
        IMPLICIT NONE
=======
                                                                 VCell,PU,PV)
      !! The subroutine compute the cofficients following Ghos Point Method        
	IMPLICIT NONE
>>>>>>> origin/CommentedCode
        TYPE(Grid),INTENT(IN):: PGrid,UGrid,VGrid
      !! The grid	
        TYPE(Cell),INTENT(IN):: PCell,UCell,VCell
      !! The cell	
        TYPE(PoissonCoefficient),INTENT(IN):: PU,PV
<<<<<<< HEAD
        REAl(dp), INTENT(in) ::Roref
=======
      !! The velocity coefficient for possion like equation 	
>>>>>>> origin/CommentedCode
        REAL(KIND=dp):: BetaP,BetaM,BetaW,BetaD,Lamda,tol
        INTEGER(kind=it4b):: i,j
      ! Set Coefficient for W,E,S,N
        tol=1.d-24
        BetaP = 1.d0/(row/Roref)
        BetaM = 1.d0/(roa/Roref)
        do i = 1,Isize
          do j = 1,Jsize
            if((PCell%vof(i,j)>=0.5d0.and.PCell%vofS(i,j)<epsi).or.            &
               (PCell%phi(i,j)<0.d0.and.PCell%vofS(i,j)>=epsi)) then ! this cell is in water and it is assigned wet cell
            ! For South Cell
              if(j>1) then
                if((PCell%vof(i,j-1)<0.5d0.and.PCell%vofS(i,j-1)<epsi).or.     &
                   (PCell%phi(i,j-1)>=0.d0.and.PCell%vofS(i,j-1)>=epsi))then ! the South Cell is in air and assigned dry cell
                  Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+            &
                                              dabs(PCell%phi(i,j-1))+tol)
                  BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP
                  PoCoef(i,j,1)=PV%Dp(i,j-1)*PGrid%dx(i,j)/VGrid%dy(i,j-1)*    &
                              BetaP*BetaM/BetaW
                elseif((PCell%vof(i,j-1)>=0.5d0.and.PCell%vofS(i,j-1)<epsi).or. &
                   (PCell%phi(i,j-1)<0.d0.and.PCell%vofS(i,j-1)>=epsi))then ! the South Cell is in water and assigned wet cell
                  BetaW=BetaM
                  PoCoef(i,j,1)=PV%Dp(i,j-1)*PGrid%dx(i,j)/VGrid%dy(i,j-1)*    &
                              BetaP*BetaM/BetaW
             ! Using the new concept to calculate the pressure gradient
                end if
              else
                BetaW=BetaM
                PoCoef(i,j,1)=PV%Dp(i,j)*PGrid%dx(i,j)/VGrid%dy(i,j)*          &
                              BetaP*BetaM/BetaW
              end if
            ! For West Cell
              if(i>1) then
                if((PCell%vof(i-1,j)<0.5d0.and.PCell%vofS(i-1,j)<epsi).or.     &
                   (PCell%phi(i-1,j)>=0.d0.and.PCell%vofS(i-1,j)>=epsi)) then ! the West Cell is in air and assigned dry cell
                  Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+            &
                                              dabs(PCell%phi(i-1,j))+tol)
                  BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP
                  PoCoef(i,j,2)=PU%Dp(i-1,j)*PGrid%dy(i,j)/UGrid%dx(i-1,j)*    &
                              BetaP*BetaM/BetaW
                elseif((PCell%vof(i-1,j)>=0.5d0.and.PCell%vofS(i-1,j)<epsi).or. &
                   (PCell%phi(i-1,j)<0.d0.and.PCell%vofS(i-1,j)>=epsi))then
                  BetaW=BetaM
                  PoCoef(i,j,2)=PU%Dp(i-1,j)*PGrid%dy(i,j)/UGrid%dx(i-1,j)*    &
                              BetaP*BetaM/BetaW
                end if
              else
                BetaW=BetaM
                PoCoef(i,j,2)=PU%Dp(i,j)*PGrid%dy(i,j)/UGrid%dx(i,j)*          &
                              BetaP*BetaM/BetaW
              end if

            ! For East Cell
              if(i<Isize) then
                if((PCell%vof(i+1,j)<0.5d0.and.PCell%vofS(i+1,j)<epsi).or.     &
                   (PCell%phi(i+1,j)>=0.d0.and.PCell%vofS(i+1,j)>=epsi)) then
                  Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+            &
                                                    dabs(PCell%phi(i+1,j))+tol)
                  BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP
                  PoCoef(i,j,3)=PU%Dp(i,j)*PGrid%dy(i,j)/UGrid%dx(i,j)*        &
                              BetaP*BetaM/BetaW
                elseif((PCell%vof(i+1,j)>=0.5d0.and.PCell%vofS(i+1,j)<epsi).or. &
                   (PCell%phi(i+1,j)<0.d0.and.PCell%vofS(i+1,j)>=epsi))then
                  BetaW=BetaM
                  PoCoef(i,j,3)=PU%Dp(i,j)*PGrid%dy(i,j)/UGrid%dx(i,j)*        &
                              BetaP*BetaM/BetaW
                end if
              else
                BetaW=BetaM
                PoCoef(i,j,3)=PU%Dp(i,j)*PGrid%dy(i,j)/UGrid%dx(i,j)*          &
                              BetaP*BetaM/BetaW
              end if
            ! For North Cell
              if(j<Jsize) then
                if((PCell%vof(i,j+1)<0.5d0.and.PCell%vofS(i,j+1)<epsi).or.     &
                   (PCell%phi(i,j+1)>=0.d0.and.PCell%vofS(i,j+1)>=epsi))  then
                  Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+            &
                                              dabs(PCell%phi(i,j+1))+tol)
                  BetaW=Lamda*BetaM+(1.d0-Lamda)*BetaP
                  PoCoef(i,j,4)=PV%Dp(i,j)*PGrid%dx(i,j)/VGrid%dy(i,j)*        &
                              BetaP*BetaM/BetaW
                elseif((PCell%vof(i,j+1)>=0.5d0.and.PCell%vofS(i,j+1)<epsi).or. &
                   (PCell%phi(i,j+1)<0.d0.and.PCell%vofS(i,j+1)>=epsi))then
                  BetaW=BetaM
                  PoCoef(i,j,4)=PV%Dp(i,j)*PGrid%dx(i,j)/VGrid%dy(i,j)*        &
                              BetaP*BetaM/BetaW
                end if
              else
                BetaW=BetaM
                PoCoef(i,j,4)=PV%Dp(i,j)*PGrid%dx(i,j)/VGrid%dy(i,j)*          &
                              BetaP*BetaM/BetaW
              end if
            elseif((PCell%vof(i,j)<0.5d0.and.PCell%vofS(i,j)<epsi).or.        &
               (PCell%phi(i,j)>=0.d0.and.PCell%vofS(i,j)>=epsi)) then ! this cell is in air and it is assigned dry cell
              if(j>1) then
              ! the South Cell is in water and assigned wet cell
                if((PCell%vof(i,j-1)>=0.5d0.and.PCell%vofS(i,j-1)<epsi).or.    &
                   (PCell%phi(i,j-1)<0.d0.and.PCell%vofS(i,j-1)>=epsi))then
                  Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+            &
                                              dabs(PCell%phi(i,j-1))+tol)
                  BetaD=Lamda*BetaP+(1.d0-Lamda)*BetaM
                  PoCoef(i,j,1)=PV%Dp(i,j-1)*PGrid%dx(i,j)/VGrid%dy(i,j-1)*    &
                              BetaP*BetaM/BetaD
                elseif((PCell%vof(i,j-1)<0.5d0.and.PCell%vofS(i,j-1)<epsi).or. &
                       (PCell%phi(i,j-1)>=0.d0.and.PCell%vofS(i,j-1)>=epsi))then ! the South Cell is in water and assigned wet cell
                  BetaW=BetaP
                  PoCoef(i,j,1)=PV%Dp(i,j-1)*PGrid%dx(i,j)/VGrid%dy(i,j-1)*    &
                              BetaP*BetaM/BetaD
                end if
              else
                BetaD=BetaP
                PoCoef(i,j,1)=PV%Dp(i,j)*PGrid%dx(i,j)/VGrid%dy(i,j)*          &
                              BetaP*BetaM/BetaD
              end if
              if(i>1) then
            ! the West Cell is in water and assigned wet cell
                if((PCell%vof(i-1,j)>=0.5d0.and.PCell%vofS(i-1,j)<epsi).or.    &
                   (PCell%phi(i-1,j)<0.d0.and.PCell%vofS(i-1,j)>=epsi))then
                  Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+            &
                                              dabs(PCell%phi(i-1,j))+tol)
                  BetaD=Lamda*BetaP+(1.d0-Lamda)*BetaM
                  PoCoef(i,j,2)=PU%Dp(i-1,j)*PGrid%dy(i,j)/UGrid%dx(i-1,j)*    &
                              BetaP*BetaM/BetaD
                elseif((PCell%vof(i-1,j)<0.5d0.and.PCell%vofS(i-1,j)<epsi).or. &
                       (PCell%phi(i-1,j)>=0.d0.and.PCell%vofS(i-1,j)>=epsi))then
                  BetaD=BetaP
                  PoCoef(i,j,2)=PU%Dp(i-1,j)*PGrid%dy(i,j)/UGrid%dx(i-1,j)*    &
                              BetaP*BetaM/BetaD
                end if
              else
                BetaD=BetaP
                PoCoef(i,j,2)=PU%Dp(i,j)*PGrid%dy(i,j)/UGrid%dx(i,j)*          &
                              BetaP*BetaM/BetaD
              end if

            ! The East Cell is in the water and assigned wet cell
              if(i<Isize) then
                if((PCell%vof(i+1,j)>=0.5d0.and.PCell%vofS(i+1,j)<epsi).or.    &
                   (PCell%phi(i+1,j)<0.d0.and.PCell%vofS(i+1,j)>=epsi))  then
                  Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+            &
                                              dabs(PCell%phi(i+1,j))+tol)
                  BetaD=Lamda*BetaP+(1.d0-Lamda)*BetaM
                  PoCoef(i,j,3)=PU%Dp(i,j)*PGrid%dy(i,j)/UGrid%dx(i,j)*        &
                              BetaP*BetaM/BetaD
                elseif((PCell%vof(i+1,j)<0.5d0.and.PCell%vofS(i+1,j)<epsi).or. &
                       (PCell%phi(i+1,j)>=0.d0.and.PCell%vofS(i+1,j)>=epsi))then
                  BetaD=BetaP
                  PoCoef(i,j,3)=PU%Dp(i,j)*PGrid%dy(i,j)/UGrid%dx(i,j)*        &
                              BetaP*BetaM/BetaD
                end if
              else
                BetaD=BetaP
                PoCoef(i,j,3)=PU%Dp(i,j)*PGrid%dy(i,j)/UGrid%dx(i,j)*          &
                              BetaP*BetaM/BetaD
              end if
           !  For North Cell
              if(j<Jsize) then
                if((PCell%vof(i,j+1)>=0.5d0.and.PCell%vofS(i,j+1)<epsi).or.    &
                   (PCell%phi(i,j+1)<0.d0.and.PCell%vofS(i,j+1)>=epsi)) then  ! the North Cell is in the Water and assigned wet cell
                  Lamda=dabs(PCell%phi(i,j))/(dabs(PCell%phi(i,j))+            &
                                              dabs(PCell%phi(i,j+1))+tol)
                  BetaD=Lamda*BetaP+(1.d0-Lamda)*BetaM
                  PoCoef(i,j,4)=PV%Dp(i,j)*PGrid%dx(i,j)/VGrid%dy(i,j)*        &
                              BetaP*BetaM/BetaD
                elseif((PCell%vof(i,j+1)<0.5d0.and.PCell%vofS(i,j+1)<epsi).or. &
                       (PCell%phi(i,j+1)>=0.d0.and.PCell%vofS(i,j+1)>=epsi))then
                  BetaD=BetaP
                  PoCoef(i,j,4)=PV%Dp(i,j)*PGrid%dx(i,j)/VGrid%dy(i,j)*        &
                              BetaP*BetaM/BetaD

                end if
              else
                BetaD=BetaP
                PoCoef(i,j,4)=PV%Dp(i,j)*PGrid%dx(i,j)/VGrid%dy(i,j)*          &
                              BetaP*BetaM/BetaD
              end if
            end if
          end do
        end do
    End Subroutine ComputePossionMatrixCoefficient

    Subroutine DeltaPressureBoundaryCondition(Proj,WB,EB,SB,NB)
	!! The subroutine insert boundary condtion for pressure
        IMPLICIT NONE
        TYPE(Projection),INTENT(INOUT):: Proj
	!! Computed pressure
        INTEGER(kind=it4b),INTENT(IN):: WB,EB,SB,NB
        !! The flag for boundary condition 
        INTEGER(kind=it4b):: i,j
        Do i = ibeg,ibeg+Isize-1
          Proj%Pp(i,jbeg-1) = dble(SB)*Proj%Pp(i,jbeg)
          Proj%Pp(i,jbeg+Jsize) = dble(NB)*Proj%Pp(i,jbeg+Jsize-1)
        End do
        Do j = jbeg,jbeg+Jsize-1
          Proj%Pp(ibeg-1,j) = dble(WB)*Proj%Pp(ibeg,j)
          Proj%Pp(ibeg+Isize,j) = dble(EB)*Proj%Pp(ibeg+Isize-1,j)
        End do
    End subroutine DeltaPressureBoundaryCondition
End module ProjectionP
