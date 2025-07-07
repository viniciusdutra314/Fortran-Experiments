      PROGRAM IntegralBoole
        IMPLICIT REAL*8(A-H,O-Z)
        a=0
        b=1
        DO i=0,9
          N=12*(2**i)
          h=(b-a)/N
          !vamos usar a regra de boole composta
          rint=7*f(a)+f(b)
          DO j=1,N,2
            rint=rint+32*f(a+h*j)
          END DO  
          DO j=2,N,4
            rint=rint+12*f(a+h*j)
          END DO
          DO j=4,N,4
            rint=rint+14*f(a+h*j)
          END DO
          rint=2*h*rint/45
          WRITE(*,*) rint       
        END DO 
      END PROGRAM


      REAL*8 FUNCTION f(x)
        REAL*8 x
        pi = 4.0d0*ATAN(1.0d0)
        f = EXP(-x) * COS(2.0d0*pi*x)
      END FUNCTION
