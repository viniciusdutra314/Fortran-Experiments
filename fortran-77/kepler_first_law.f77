      PROGRAM PrimeiraLeiDeKepler
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(UA=1.496d11)
      character name*100
      character dir_position*100
      character dir_distances*100
      DIMENSION p(2)
      DIMENSION v(2)
      WRITE(*,*) 'Primeira lei de Kepler (em SI)'
      WRITE(*,*) 'Digite x0,y0,v0,vy,tf,dt,nome_arquivo'
      READ(*,*) p(1),p(2),v(1),v(2),tf,dt,name
      WRITE(dir_position,30) name
  30  FORMAT("planetas_posicoes/",(A10),'.csv')
      WRITE(dir_distances,40) name
  40  FORMAT("planetas_distancias/",(A10),'.csv')
      t=0
      OPEN(unit=10,file=dir_position)
      d_max=SQRT(p(1)**2+p(2)**2)
      d_min=SQRT(p(1)**2 + p(2)**2)
      x_perielio=0
      x_afelio=0
      y_perielio=0
      y_afelio=0
      DO WHILE(t<tf)
        t=t+dt
        a=gravity(p(1),p(2))
        p(1)=p(1)+v(1)*dt+(a*p(1)/2)*dt**2
        p(2)=p(2)+v(2)*dt+(a*p(2)/2)*dt**2
        v(1)=v(1)+a*p(1)*dt
        v(2)=v(2)+a*p(2)*dt
        WRITE(10,*) p(1), ',', p(2)
        distancia=SQRT(p(1)**2 + p(2)**2)
        IF (d_max < distancia) THEN
          d_max=distancia
          x_afelio=p(1)
          y_afelio=p(2)
        END IF
        IF (d_min > distancia) THEN
          d_min=distancia
          x_perielio=p(1)
          y_perielio=p(2)
        END IF
      END DO
      a=(d_max+d_min)/2
      c=a-d_min
      CLOSE(10)
      OPEN(unit=10,file=dir_position)
      OPEN(unit=20,file=dir_distances)
      t=0
      IF (x_perielio.ne.0) THEN
        theta=ATAN(y_perielio/x_perielio)
      ELSE
        theta=ATAN(y_afelio,x_afelio)
      END IF
      rx=2*c*COS(theta)
      ry=2*c*SIN(theta)
      DO WHILE(t<tf)
        t=t+dt
        READ(10,*) x,y
        distance=SQRT((x)**2+(y)**2)
        distance=distance+SQRT((x-rx)**2+(y-ry)**2)
        WRITE(20,*) distance
      END DO
      CLOSE(10)
      CLOSE(20)
      END PROGRAM

      REAL*8 FUNCTION gravity(x,y)
        IMPLICIT REAL*8(A-H,O-Z)
        PARAMETER(G=6.67458d-11,rm_sol=1.989d30)
        gravity=-G*rm_sol/(sqrt(x**2+y**2))**3 
      END FUNCTION