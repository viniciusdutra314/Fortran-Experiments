module rotinas
   use iso_fortran_env
   implicit none
   abstract interface
   function beta_function(time_step) result(beta)
      import :: real64, int64
      integer(int64), intent(in) :: time_step
      real(real64) :: beta
      end function beta_function
   end interface

contains
   subroutine random_initial_condition(L,spins)
      integer(int64), intent(in) :: L
      integer(int8),intent(inout) :: spins(L,L)
      real :: random_numbers(L,L)
      call random_number(random_numbers)
      where (random_numbers < 0.5) spins = -1
      where (random_numbers >= 0.5) spins = 1
   end subroutine

   function bernoulli_trial(probability) result(number)
      real(real64), intent(in) :: probability
      logical :: number
      real :: rand
      call random_number(rand)
      if (rand < probability) then
         number = .true.
      else
         number = .false.
      end if
   end function bernoulli_trial

   subroutine random_point(L, x_cor, y_cor)
      integer(int64), intent(in) :: L
      integer(int64), intent(out) :: x_cor, y_cor
      real(real64) :: rand1, rand2
      
      call random_number(rand1)
      call random_number(rand2)
      x_cor = NINT(rand1 * (L-1)) + 1
      y_cor = NINT(rand2 * (L-1)) + 1
   end subroutine random_point

   function total_energy(spins,iplus,iminus) result(energy)
      integer(int8), intent(in) :: spins(:,:)
      integer(int64), intent(in) :: iplus(:), iminus(:)
      integer(int64) :: energy,i,j
      energy = 0
      do i = 1, size(spins,1)
         do j = 1, size(spins,2)
            energy = energy - spins(i,j) * (spins(iminus(i),j)+spins(iplus(i),j)+spins(i,iminus(j))+spins(i,iplus(j)))
         end do
      end do
      energy = energy / 2
   end function total_energy

   subroutine probabilities_table(probabilities, beta)
      real(real64), intent(inout) :: probabilities(-4:4)
      real(real64), intent(in) :: beta
      integer(int64) :: i
      do i = -4,4
             probabilities(i) = exp(-beta*i)
      end do   
   end subroutine probabilities_table


   subroutine simulation(L,beta_func,filename,initial_condition,mc_steps,write_frequency)
      procedure(beta_function) :: beta_func
      integer(int64), intent(in) :: L,mc_steps
      character(len=*), intent(in) :: filename,initial_condition,write_frequency
      
      integer(int8) spins(L,L)
      integer(int64) :: iplus(L),iminus(L),i,j,x_cor,y_cor,delta_E,E,M,N
      real(real64) beta,probabilities(-4:4)
      !criando condições de contorno periódicas
      iplus = [(i+1, i=1,L)]
      iplus(L) = 1
      iminus = [(i-1, i=1,L)]
      iminus(1) = L
      
      select case (trim(initial_condition))
      case ("random")
         call random_initial_condition(L, spins)
      case ("all_up")
         spins = 1
      case ("all_down")
         spins = -1
      case default
         spins = 1  
      end select
      
      E= total_energy(spins,iplus,iminus)
      M=0 ! Poderiamos usar sum(spins), mas isso causaria um overflow 
      DO i = 1, L
         do j = 1, L
            M = M + spins(i,j)
         end do
      end do

      N = size(spins)
      open(10,file="raw_data/" // filename,status='replace')
      DO i = 1 , mc_steps
         beta = beta_func(i)
         call probabilities_table(probabilities, beta)
         if (write_frequency=="outer") then
            write(10,*) M/real(N),E/real(N),beta
         endif

         DO j = 1 , 10*size(spins) 
            call random_point(L, x_cor, y_cor)
            delta_E= spins(x_cor,y_cor)*(spins(iminus(x_cor),y_cor)+spins(iplus(x_cor),y_cor)+&
            spins(x_cor,iminus(y_cor))+spins(x_cor,iplus(y_cor)))

            if (delta_E<0 .or. bernoulli_trial(probabilities(delta_E))) then
               spins(x_cor,y_cor) = - spins(x_cor,y_cor)
               E=E + 2*delta_E 
               M= M + 2 * spins(x_cor,y_cor)
            end if

            if (write_frequency=="inner") then
                write(10,*) M/real(N),E/real(N),beta
            end if
         END DO
      END DO
      close(10)

      open(20,file="raw_data/final_"//trim(filename),status='replace')
      write(20,*) spins
      close(20)

   end subroutine simulation

   function beta_pequeno(time_step) result(beta)
      integer(int64), intent(in) :: time_step
      real(real64) :: beta
      beta = 0.1
   end function beta_pequeno

   function beta_3(time_step) result(beta)
      integer(int64), intent(in) :: time_step
      real(real64) :: beta
      beta = 3
   end function beta_3

   function beta_meio(time_step) result(beta)
      integer(int64), intent(in) :: time_step
      real(real64) :: beta
      beta = 0.5
   end function beta_meio

   function beta_recozimento(time_step) result(beta)
      integer(int64), intent(in) :: time_step
      real(real64) :: beta
      beta = 0 + 0.001*time_step
      if (beta > 3) then
         beta = 3
      end if
   end function beta_recozimento


   function beta_histerese_curta(time_step) result(beta)
      integer(int64), intent(in) :: time_step
      real(real64) :: beta
      beta = 0.001*time_step
      if (beta > 1.75) then
         beta = 1.75 - 0.001*(time_step-1.75/0.001) 
      end if
   end function beta_histerese_curta

    function beta_histerese_longa(time_step) result(beta)
      integer(int64), intent(in) :: time_step
      real(real64) :: beta
      beta = 0.0001d0*time_step
      if (beta > 1.75) then
         beta = 1.75 - 0.0001d0*(time_step-1.75/0.0001d0) 
      end if
   end function beta_histerese_longa


end module rotinas

program simulation_monte_carlo
   use iso_fortran_env   
   use rotinas
   call simulation(100_int64, beta_pequeno, "magnetization_beta_01_L100.data", "all_up",20_int64,"inner")
   call simulation(60_int64, beta_pequeno, "magnetization_beta_01_L60.data", "all_up",20_int64,"inner")
   call simulation(60_int64, beta_3, "magnetization_beta_3_L60.data", "all_up",20_int64,"inner")
   call simulation(100_int64, beta_3, "magnetization_beta_3_L100.data", "all_up",20_int64,"inner")
   call simulation(100_int64, beta_recozimento, "energy_recozimento.data", "random",3000_int64,"outer")
   call simulation(100_int64, beta_3, "energy_tempera.data", "random",3000_int64,"outer")
   call simulation(60_int64, beta_histerese_curta, "L_60_histerese_0001.data", "random",3500_int64,"outer")
   call simulation(60_int64, beta_histerese_longa, "L_60_histerese_00001.data", "random",35000_int64,"outer")
   call simulation(80_int64, beta_histerese_curta, "L_80_histerese_0001.data", "random",3500_int64,"outer")
   call simulation(80_int64, beta_histerese_longa, "L_80_histerese_00001.data", "random",35000_int64,"outer")
   call simulation(100_int64, beta_histerese_curta, "L_100_histerese_0001.data", "random",3500_int64,"outer")
   call simulation(100_int64, beta_histerese_longa, "L_100_histerese_00001.data", "random",35000_int64,"outer")
   
   call simulation(4_int64, beta_meio, "quebra_simetria_4.data", "random",2000_int64,"outer")
   call simulation(5_int64, beta_meio, "quebra_simetria_5.data", "random",2000_int64,"outer")
   call simulation(6_int64, beta_meio, "quebra_simetria_6.data", "random",2000_int64,"outer")
   call simulation(7_int64, beta_meio, "quebra_simetria_7.data", "random",2000_int64,"outer")
   call simulation(8_int64, beta_meio, "quebra_simetria_8.data", "random",2000_int64,"outer")
   call simulation(9_int64, beta_meio, "quebra_simetria_9.data", "random",2000_int64,"outer")
   call simulation(10_int64, beta_meio, "quebra_simetria_10.data", "random",2000_int64,"outer")

end program simulation_monte_carlo