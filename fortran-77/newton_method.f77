      PROGRAM newton
      IMPLICIT REAL*8(A-H,O-Z)
      epsilon=1e-6
      WRITE(*,*) 'Escolha um valor inicial x0'
      READ(*,*) x
  100 FORMAT(A12,A12,A12)
  200 FORMAT(I12,1PE12.3,1PE12.3)
      WRITE(*,100) 'Iteração','x', 'p(x)'
      WRITE(*,200) 0,x,p(x)
      iteracao=0
      DO WHILE(ABS(rnewton_step).gt.epsilon.or.iteracao.eq.0 )
        rnewton_step=p(x)/diff(x)
        x=x-rnewton_step
        iteracao=iteracao+1
        WRITE(*,200) iteracao, x,p(x)
      END DO
      END PROGRAM


      REAL*8 FUNCTION p(x)
        REAL*8 x
        p=x**3 -4*x**2 -59*x +126
      END FUNCTION

      REAL*8 FUNCTION diff(x)
        REAL*8 x
        diff=3*x**2 -8*x -59
      END FUNCTION
