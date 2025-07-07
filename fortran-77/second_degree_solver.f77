      PROGRAM Bhaskara
      WRITE (*,*) "Insira os coeficientes a,b,c:"
      READ * , a,b,c
      delta = b*b -4*a*c
      IF (delta>0) THEN
        WRITE (*,*) "Existem duas raizes reais"
        araiz=(-b +SQRT(delta))/(2*a)
        braiz=(-b -SQRT(delta))/(2*a)
        WRITE(*,*) araiz, braiz
      ELSE IF (delta<0) THEN
        WRITE (*,*) "Não existem raizes reais"
      ELSE
        WRITE (*,*) "As raizes são repetidas"
        WRITE (*,*)  -b/(2*a)
      END IF
      END PROGRAM
