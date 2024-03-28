! Dynamical instabilities cause extreme events in a theoretical Brusselator model.
! Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
! Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)

IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION X(10), Y(10), AK(10,4), F(10)
COMMON B

! Open output files for writing
open(8, file='bo_phase.dat')
open(9, file='bo_poin.dat')

! Initial conditions
T = 0.0D0
X(1) = 0.01
X(2) = 0.02

NH = 1098
W1 = 0.7D0
PI = 4.0 * ATAN(1.0)
H = (2.0 * PI) / (W1 * FLOAT(NH))

! Loop for simulating the dynamical system
DO J = 1, 70000000
    T = T + H
    
    ! Call RK4 subroutine to perform integration
    CALL RK4(X, T, H, F)
    
    ! Write data to output file 9 ('bo_poin.dat')
    IF (J > 32000000 .AND. F(2) < 0.0D0 .AND. F1 > 0.0D0) THEN
        WRITE(9, *) T, X(1), X(2), AM * DSIN(W1 * T)
    END IF
    
    ! Update F1 and F2 for the next iteration
    F2 = F1
    F1 = F(1)
END DO

STOP
END

! Subroutine to perform the fourth-order Runge-Kutta method
SUBROUTINE RK4(X, T, H, F)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION X(10), Y(10), AK(10,4), F(10)
N1 = 2

DO I = 1, N1
    Y(I) = X(I)
END DO

CALL SUB(T, Y, F)
DO I = 1, N1
    AK(I,1) = H * F(I)
    Y(I) = X(I) + AK(I,1) / 2.0D0
END DO

CALL SUB(T + H / 2.0D0, Y, F)
DO I = 1, N1
    AK(I,2) = H * F(I)
    Y(I) = X(I) + AK(I,2) / 2.0D0
END DO

CALL SUB(T + H / 2.0D0, Y, F)
DO I = 1, N1
    AK(I,3) = H * F(I)
    Y(I) = X(I) + AK(I,3)
END DO

CALL SUB(T + H, Y, F)
DO I = 1, N1
    AK(I,4) = H * F(I)
    X(I) = X(I) + 1.0D0 / 6.0D0 * (AK(I,1) + 2.0D0 * (AK(I,2) + AK(I,3)) + AK(I,4))
END DO

RETURN
END

! Subroutine to define the equations of the system
SUBROUTINE SUB(T, X, F)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION X(10), F(10)
COMMON B

A = 0.2
AM = 0.06
W1 = 0.7

F(1) = A - B * X(1) - X(1) + X(1) * X(1) * X(2) + AM * DSIN(W1 * T)
F(2) = B * X(1) - X(1) * X(1) * X(2)

RETURN
END
