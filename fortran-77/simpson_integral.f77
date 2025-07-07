      PROGRAM IntegralSimpson
      IMPLICIT REAL*8(A-H,O-Z)
      a=0
      b=1
      DO i=0,9
        N=12*(2**i)
        h=(b-a)/N
        rint=0
        DO j=1,N-1
            x_a=a+h*(j-1)
            x_b=a+h*j
            rint=rint+((x_b-x_a)/6)*(f(x_a)+4*f((x_a+x_b)/2)+f(x_b))
        END DO
        WRITE(*,*) rint        
        END DO 
        END PROGRAM


      REAL*8 FUNCTION f(x)
        REAL*8 x
        pi = 4.0d0*ATAN(1.0d0)
        f = EXP(-x) * COS(2.0d0*pi*x)
      END FUNCTION
