      PROGRAM lei_dos_periodos
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(UA=1.496d11)
      character name*50
      DIMENSION p(2)
      DIMENSION v(2)
      WRITE(*,*) 'Terceira Lei de Kepler (dados em SI)'
      WRITE(*,*) 'Digite x0,y0,v0,vy,tf,dt'
      READ(*,*) p(1),p(2),v(1),v(2),tf,dt
      t=0
      sign=1
      sign_changes=0
      period=0
      d_max=SQRT(P(1)**2 + P(2)**2)
      d_min=SQRT(P(1)**2 + P(2)**2)
      DO WHILE(t<tf)
        t=t+dt
        a=gravity(p(1),p(2))
        p(1)=p(1)+v(1)*dt+(a*p(1)/2)*dt**2
        p(2)=p(2)+v(2)*dt+(a*p(2)/2)*dt**2
        v(1)=v(1)+a*p(1)*dt
        v(2)=v(2)+a*p(2)*dt
        distance=SQRT(P(1)**2 + P(2)**2)
        IF (distance>d_max) THEN
          d_max=distance
        END IF
        IF (distance<d_min) THEN
          d_min=distance
        END IF
        IF ((p(1)/ABS(p(1))).ne.sign) THEN
          sign_changes=sign_changes+1
          IF (sign_changes.eq.2) THEN
            a=(d_max+d_min)/2
            WRITE(*,*) 'Semieixo-maior e periodo'
            WRITE(*,*)  a, 2*t
            t=tf+dt
          END IF
        END IF
      END DO
      END PROGRAM

      REAL*8 FUNCTION gravity(x,y)
        IMPLICIT REAL*8(A-H,O-Z)
        G=6.67458d-11
        rm_sol=1.989d30
        gravity=-G*rm_sol/(sqrt(x**2+y**2))**3 
      END FUNCTION