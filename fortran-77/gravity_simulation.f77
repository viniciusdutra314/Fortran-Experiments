      PROGRAM mudanca_referencial
      IMPLICIT REAL*8(A-H,O-Z)
      character name_a*100
      character name_b*100
      character dir*100
      DIMENSION p(2,2)
      DIMENSION v(2,2)
      WRITE(*,*) 'Dados do primeiro planeta:'
      WRITE(*,*) 'Digite x0,y0,vx,vy,nome'
      READ(*,*) p(1,1),p(1,2),v(1,1),v(1,2),name_a
      WRITE(*,*) 'Dados do segundo planeta:'
      READ(*,*) p(2,1),p(2,2),v(2,1),v(1,2),name_b
      dir='dados//' // TRIM(name_a) // '-' // TRIM(name_b) // '.csv'
      WRITE(*,*) 'Dados da simulação tf,dt'
      READ(*,* ) tf,dt
      t=0
      OPEN(unit=10,file=dir)
      WRITE(10,100) 
  100 FORMAT('x_values,y_values,distance')
      DO WHILE(t<tf)
        t=t+dt
        a1=gravity(p(1,1),p(1,2))
        p(1,1)=p(1,1)+v(1,1)*dt+(a1*p(1,1)/2)*dt**2
        p(1,2)=p(1,2)+v(1,2)*dt+(a1*p(1,2)/2)*dt**2
        v(1,1)=v(1,1)+a1*p(1,1)*dt
        v(1,2)=v(1,2)+a1*p(1,2)*dt
        a2=gravity(p(2,1),p(2,2))
        p(2,1)=p(2,1)+v(2,1)*dt+(a2*p(2,1)/2)*dt**2
        p(2,2)=p(2,2)+v(2,2)*dt+(a2*p(2,2)/2)*dt**2
        v(2,1)=v(2,1)+a2*p(2,1)*dt
        v(2,2)=v(2,2)+a2*p(2,2)*dt
        diff_x=p(2,1)-p(1,1)
        diff_y=p(2,2)-p(1,2)
        distance=SQRT(diff_x**2 + diff_y**2)
        WRITE(10,*) diff_x, ',',diff_y,',', distance
      END DO
      CLOSE(10)
      END PROGRAM

      REAL*8 FUNCTION gravity(x,y)
        IMPLICIT REAL*8(A-H,O-Z)
        PARAMETER(G=6.67458d-11,rm_sol=1.989d30)
        gravity=-G*rm_sol/(sqrt(x**2+y**2))**3 
      END FUNCTION