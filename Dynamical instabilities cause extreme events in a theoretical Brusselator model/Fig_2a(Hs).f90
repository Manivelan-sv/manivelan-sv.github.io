
! Dynamical instabilities cause extreme events in a theoretical Brusselator model.
! Manivelan, S. V., Sabarathinam, S., Thamilmaran, K., & Manimehan, I.
! Chaos, Solitons & Fractals, 180, 114582. (https://doi.org/10.1016/j.chaos.2024.114582)


	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(10),Y(10),AK(10,4),F(10),Xn(10),Xm(10),ab(10000000),SX(25000000)
	COMMON B
	INTEGER :: tt
!	*******Output files*********
	open(8,file='hs.dat')
!	*******Initial conditions***********
	T=0.0D0
   
!	*******Step size************
     	A = 0.219D0
		B = 0.9D0
	 
110     B = B+0.00001
      	
	X(1)=0.01D0
	X(2)=0.02D0 

	AM=0.06
	W1=0.7 
	PI = 4.0*ATAN(1.0)	
	H=(2.0*PI)/(W1*FLOAT(1098))

	nt = 4
    
    SS = 0.0D0
    tt = 0
    SD = 0.0D0
    XAVG = 0.0D0
	XMAX = -1.0D0

	
!       *********Phase loop********2*10^5***
	
	DO J=1,2000000
		T=T+H
		CALL RK4(X,T,H,F)
		f2=f1
		f1=f(1)
		IF(J.GT. 20000)then
		if(f1.lt.0.0.and.f2.gt.0.0) then
			si = abs(X(1))
			SS = SS + si
                	tt = tt + 1
			SX(tt) = si
                	IF (si > XMAX) THEN
                    		XMAX = si
                	END IF
		END IF
		END IF
	END DO 
	
!!!!!!!!!!!!!!!!!!!! Standard Deviation Calculation !!!!!!!!!!!!!!!!!!!!!!!!
	XAVG = SS / tt
    	CUM = 0.0D0
   	DO J = 1, tt
        	CUM = CUM + (SX(J) - XAVG)**2
    	END DO
    	SD = SQRT(CUM / (tt))
    
    Hs = XAVG + (nt * SD)
	Hs5 = XAVG + (5 * SD)
	Hs6 = XAVG + (6 * SD)
	Hs7 = XAVG + (7 * SD)
	Hs8 = XAVG + (8 * SD)
	SEE = XAVG + (12.0D0 * SD)
	d=(XMAX - XAVG)/SD
	WRITE(8,*)B,Hs,Hs5,Hs6,Hs7,Hs8,SEE,d

	IF(B.lT. 1.2D0) then
		GOTO 110
	endif
    STOP
	END
!	******************* RK4 subroutine ********
	SUBROUTINE RK4(X,T,H,F)	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(10),Y(10),AK(10,4),F(10)
!
	N1=2 
       
! 
	DO I=1,N1
		Y(I) = X(I)
	END DO
!	
	CALL SUB(T,Y,F)
	DO I=1,N1
	AK(I,1)=H*F(I)
	Y(I)=X(I)+AK(I,1)/2.0D0	
	END DO	
!	
	CALL SUB(T+H/2.0,Y,F)

	DO I=1,N1
		AK(I,2)=H*F(I)
		Y(I)=X(I)+AK(I,2)/2.0D0
	END DO
!C
	CALL SUB(T+H/2.0,Y,F)
	DO I=1,N1
		AK(I,3)=H*F(I)
		Y(I)=X(I)+AK(I,3)
	END DO	
!C
	CALL SUB(T+H,Y,F)
	DO I=1,N1
		AK(I,4)=H*F(I)
		X(I)=X(I)+1/6.0D0*(AK(I,1)+2.0D0*(AK(I,2)+AK(I,3))+AK(I,4))
	END DO

!C                 
	RETURN
	END	                                                                                         
!c
!C       ***************Subroutine for Equations ********************
!C
	SUBROUTINE SUB(T,X,F)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION X(10),F(10)
	COMMON B

	AM=0.06
	W1=0.7

        F(1) = A-B*X(1)-X(1)+X(1)*X(1)*X(2)+AM*DSIN(W1*T)
       	F(2) = B*X(1)-X(1)*X(1)*X(2)

	RETURN
	END
