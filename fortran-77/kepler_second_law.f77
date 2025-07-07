      PROGRAM SegundaLeiDeKepler
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(Ua=1.496d11)
      character name*50
      character diretorio*50
      DIMENSION p(2)
      DIMENSION v(2)
      WRITE(*,*) 'Segunda lei de Kepler (dados em SI)'
      WRITE(*,*) 'Digite x0,y0,v0,vy,tf,dt,nome_arquivo'
      READ(*,*) p(1),p(2),v(1),v(2),tf,dt,name
      WRITE(diretorio,30) name
  30  FORMAT("planetas_dados/",(A10),'.csv')
      t=0
      OPEN(unit=10,file=diretorio)
      WRITE(10,100) 
  100 FORMAT('momento_angular')
      DO WHILE(t<tf)
        t=t+dt
        a=gravity(p(1),p(2))
        p(1)=p(1)+v(1)*dt+(a*p(1)/2)*dt**2
        p(2)=p(2)+v(2)*dt+(a*p(2)/2)*dt**2
        v(1)=v(1)+a*p(1)*dt
        v(2)=v(2)+a*p(2)*dt
        area_time=ABS(p(1)*v(2)-v(1)*p(2))/2
        WRITE(10,*) area_time
      END DO
      WRITE(*,*) 'Dados salvos em ' , diretorio
      CLOSE(10)
      END PROGRAM
      REAL*8 FUNCTION gravity(x,y)
        IMPLICIT REAL*8(A-H,O-Z)
        PARAMETER(G=6.67458d-11,rm_sol=1.989d30)

        gravity=-G*rm_sol/(sqrt(x**2+y**2))**3 
      END FUNCTION