subroutine gaussiana(array_size,y,mu,sigma,dx)
   integer, intent(in) :: array_size
   real(kind=8), intent(out) :: y(array_size)
   real(kind=8), intent(in) :: mu,sigma,dx
   do i=1,array_size
      y(i)=exp(-((i*dx-mu)/sigma)**2)
   end do
end subroutine

subroutine propagate_and_save(filename_output,dx,L,c,idownsampling_factor,r,time_simulated)
   character(len=*) filename_output
   real(kind=8),intent(in) :: dx,r,c,L,time_simulated
   integer,intent(in):: idownsampling_factor
   real(kind=8), allocatable :: y_futuro(:),y_presente(:),y_passado(:),temp(:)
   Nt=int(time_simulated/(r*dx/c))
   Nx=int(L/dx)
   allocate(y_futuro(Nx))
   allocate(y_presente(Nx))
   allocate(y_passado(Nx))
   allocate(temp(Nx))
   !condicoes inicias
   call gaussiana(Nx,y_passado,L/3.d0,L/30.d0,dx)

   y_presente=y_passado !velocidade zero
   !condicoes de contorno
   y_futuro(1)=0;y_futuro(Nx)=0
   y_presente(1)=0;y_presente(Nx)=0
   y_passado(1)=0;y_passado(Nx)=0

   open(unit=10, file=filename_output,status="unknown",action="write")
   write(10,*) y_passado
   write(10,*) y_presente
   DO it=1,Nt
      DO ix=1+1,Nx-2
         y_futuro(ix) = 2*(1-r**2)*y_presente(ix) + r*r*(y_presente(ix+1) + y_presente(ix-1)) - y_passado(ix)
      END DO
      temp=y_passado
      y_passado=y_presente
      y_presente=y_futuro
      y_futuro=temp
      if (mod(it,idownsampling_factor)==0) then
         write (10,*) y_futuro
      end if
   END DO

   close(10)

end subroutine

program propagacao_gaussiana
   implicit real(kind=8) (a-h,o-z)
   real(kind=8) L
   dx=0.001d0
   c=300.d0
   L=1
   idownsampling_factor=10
   call propagate_and_save("data/gaussiana_r1.data",dx,L,c,idownsampling_factor,1.0d0,0.006d0)
   call propagate_and_save("data/achar_periodo_gaussiana_r1.data",dx,L,c,idownsampling_factor,1.0d0,0.03d0)
   call propagate_and_save("data/gaussiana_r2.data",dx,L,c,idownsampling_factor,2.0d0,0.006d0)
   call propagate_and_save("data/gaussiana_r0.25.data",dx,L,c,idownsampling_factor,0.25d0,0.006d0)

end program propagacao_gaussiana
