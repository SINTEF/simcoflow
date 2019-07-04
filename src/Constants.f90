Module Constants
    USE PrecisionVar
    IMPLICIT NONE
    PRIVATE
    REAL(KIND=dp),PARAMETER,PUBLIC::pi=4.d0*datan(1.d0),Cp=1.005d3,            &
                                    kT=0.0271d0,kTw=0.0271d0,g = 9.81d0,factor=0.5d0
    ! Not really constants, but they are defined as parameters for now

    REAL(KIND=dp),PUBLIC,PARAMETER::nuw=1.0034d-6,nua=1.506d-5,                &
                                    roa=1.225d0,row=998.2d0 ! At T = 20oC
    REAL(KIND=dp),PUBLIC,PARAMETER::epsi=1.d-3,epsiF=1.d-2,BetaVis=0.5d0
End Module Constants

