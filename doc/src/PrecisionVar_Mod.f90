Module PrecisionVar
    IMPLICIT NONE
    INTEGER,PARAMETER:: sp = selected_real_kind(6,35)   ! 10 digit after ","  and form 10^-35 to 10^35
    INTEGER,PARAMETER:: dp = selected_real_kind(15,307) ! 20 digit after "," and from 10^-307 to 10^307
    INTEGER,PARAMETER:: it1b = selected_int_kind(2) ! from -10^2 to 10^2
    INTEGER,PARAMETER:: it2b = selected_int_kind(3) ! from -10^3 to 10^3
    INTEGER,PARAMETER:: it4b = selected_int_kind(6) ! from -10^6 to 10^6
    INTEGER,PARAMETER:: it8b = selected_int_kind(12) ! from -10^12 to 10^12
End module
