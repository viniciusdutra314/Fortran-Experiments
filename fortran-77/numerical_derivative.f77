      PROGRAM derivada
      IMPLICIT REAL*8(A-H,O-Z)
      dimension h_list(19)
      parameter(df=9.79678201383810,ddf=64.0983245494792,
     &dddf=671.514613457867)
      h_list=[5e-1,2e-1,1e-1,5e-2,1e-2,5e-3,1e-3,5e-4,  
     &1e-4,5e-5,1e-5,1e-6,1e-7,1e-8,1e-10,1e-15,
     &1e-16,1e-17,1e-20]
      x=0.5
  100 FORMAT(1PE12.2,1PE12.2,1PE12.2,1PE12.2,1PE12.2,1PE12.2,1PE12.2)
  200 FORMAT(A12,A12,A12,A12,A12,A12,A12)    
      WRITE(*,200) 'h', 'd3s', 'd2frente','d2traz','d5s','dd5s','ddd5as'
      DO i=1,size(h_list) 
        h=h_list(i)
        d3s=(f(x+h)-f(x-h))/(2*h) -df
        d2frente=(f(x+h)-f(x))/(h) -df
        d2traz=(f(x)-f(x-h))/(h) -df
        d5s=(f(x-2*h)-8*f(x-h)+8*f(x+h)-f(x+2*h))/(12*h) -df
        dd5s=(-f(x-2*h)+16*f(x-h)-30*f(x)+16*f(x+h)-f(x+2*h))/(12*h**2)
        dd5s=dd5s -ddf
        ddd5as=(-f(x-2*h)+2*f(x-h)-2*f(x+h)+f(x+2*h))/(2*h**3) -dddf
        WRITE(*,100) h, d3s,d2frente,d2traz,d5s, dd5s ,ddd5as
      END DO
      END PROGRAM

      REAL*8 FUNCTION f(x)
        REAL*8 x
        f=EXP(x/2)*TAN(2*x)
      END FUNCTION
