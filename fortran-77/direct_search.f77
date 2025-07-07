      PROGRAM Busca
      IMPLICIT REAL*8 (A-H,O-Z)
      epsilon=1e-6
      times=0
      WRITE(*,*) 'Escolha um intervalo [a,b]'
      READ(*,*) a,b
      WRITE(*,*) 'Numéro máximo de iterações dx=(b-a)/max_iteracoes'
      READ(*,*) max_iteracoes
      dx=(b-a)/max_iteracoes
      past_func_value=p(a)
      DO WHILE(max_iteracoes>times)
        x=a+dx*times
        times=times+1
        IF(ABS(p(x)).lt.epsilon) THEN
          WRITE(*,*) 0, "Raiz encontrada", x
        END IF
        IF (p(x)*past_func_value.lt.0.and.ABS(p(x)).gt.epsilon) THEN
          iteracoes=1
          WRITE(*,*) "         ", "----Bissecção Ativada-----"
          WRITE(*,*) iteracoes, x
          a_bis=x-dx
          b_bis=x
          c=(a_bis+b_bis)/2
          DO WHILE(ABS(p(c))>epsilon)
            iteracoes=iteracoes+1
            c=(a_bis+b_bis)/2
            WRITE(*,*) iteracoes, c
            IF (p(a_bis)*p(c)>0) THEN
               a_bis=c
            ELSE 
               b_bis=c
            END IF
           END DO
           WRITE(*,*) iteracoes+1,"Raiz encontrada", c
         END IF
        past_func_value=p(x) 
      END DO
      END PROGRAM

      REAL*8 FUNCTION p(x)
        REAL*8 x
        p=x**3 -4*x**2 -59*x +126
      END FUNCTION
