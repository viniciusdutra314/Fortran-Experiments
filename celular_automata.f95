module simulation_module
   implicit none
contains
   function integer_to_binary(n) result(binary)
      integer, intent(in) :: n
      logical :: binary(8)
      integer :: i
      do i = 0, 7
         binary(i+1) = mod(n/(2**i),2)==1
      end do
   end function integer_to_binary

   function apply_rule(left, center, right, rule_number) result(res)
      logical, intent(in) :: left, center, right
      logical, intent(in) :: rule_number(8)
      logical :: res
      integer :: index

      index = 0
      if (right) index = index + 1
      if (center) index = index + 2
      if (left) index = index + 4

      res = rule_number(index+1)
   end function apply_rule


   subroutine simulate_ACD(initial_condition, rule_number, iterations, size, filename_output)
      integer, intent(in) :: iterations, size,rule_number
      logical :: rule_number_binary(8)
      logical, dimension(size), intent(in) :: initial_condition
      logical, dimension(size) :: array, copy_array
      integer :: i, ileft, iright, it
      character(len=*), intent(in) :: filename_output
      rule_number_binary =integer_to_binary(rule_number)
      array = initial_condition
      open(unit=10, file=filename_output, status="unknown", action="write")
      do it=1, iterations
         write(10, *) array
         copy_array = array
         do i = 1, size
            ileft = modulo(i-2, size) + 1
            iright = modulo(i, size) + 1
            array(i) = apply_rule(copy_array(ileft), copy_array(i), copy_array(iright), rule_number_binary)
         end do
      end do
      close(10)
   end subroutine simulate_ACD
end module simulation_module

program tarefa1
   use simulation_module
   implicit none
   integer, parameter :: L_alcaraz=100, iterations_alcaraz=20
   integer, parameter :: L_wolphram=100, iterations_wolphram=50
   logical, dimension(L_alcaraz) :: c1, c2, c3
   logical, dimension(L_wolphram) :: c4
   real :: x(L_alcaraz)

   integer, parameter :: RULE_MAIORIA = 232
   integer, parameter :: RULE_EPIDEMIA = 254
   integer, parameter :: RULE_CONTRA = 51

   c1 = .false.
   c2 = .true.
   call random_number(x)
   where (x > 0.5)
      c3 = .false.
   elsewhere
      c3 = .true.
   end where

   call simulate_ACD(c1, RULE_MAIORIA, iterations_alcaraz, L_alcaraz, "dados/c1_regra_da_maioria.data")
   call simulate_ACD(c2, RULE_MAIORIA, iterations_alcaraz, L_alcaraz, "dados/c2_regra_da_maioria.data")
   call simulate_ACD(c3, RULE_MAIORIA, iterations_alcaraz, L_alcaraz, "dados/c3_regra_da_maioria.data")

   call simulate_ACD(c1, RULE_CONTRA, iterations_alcaraz, L_alcaraz, "dados/c1_regra_do_contra.data")
   call simulate_ACD(c2, RULE_CONTRA, iterations_alcaraz, L_alcaraz, "dados/c2_regra_do_contra.data")
   call simulate_ACD(c3, RULE_CONTRA, iterations_alcaraz, L_alcaraz, "dados/c3_regra_do_contra.data")

   call simulate_ACD(c1, RULE_EPIDEMIA, iterations_alcaraz, L_alcaraz, "dados/c1_regra_da_epidemia.data")
   call simulate_ACD(c2, RULE_EPIDEMIA, iterations_alcaraz, L_alcaraz, "dados/c2_regra_da_epidemia.data")
   call simulate_ACD(c3, RULE_EPIDEMIA, iterations_alcaraz, L_alcaraz, "dados/c3_regra_da_epidemia.data")
   c4 = .false.
   c4(L_wolphram/2+1) = .true.
   call simulate_ACD(c4, 30, iterations_wolphram, L_wolphram, "dados/c4_regra_30.data")
   call simulate_ACD(c4, 54, iterations_wolphram, L_wolphram, "dados/c4_regra_54.data")
   call simulate_ACD(c4, 60, iterations_wolphram, L_wolphram, "dados/c4_regra_60.data")
   call simulate_ACD(c4, 90, iterations_wolphram, L_wolphram, "dados/c4_regra_90.data")
   call simulate_ACD(c4, 126, iterations_wolphram, L_wolphram, "dados/c4_regra_126.data")
   call simulate_ACD(c4, 150, iterations_wolphram, L_wolphram, "dados/c4_regra_150.data")
   call simulate_ACD(c4, 188, iterations_wolphram, L_wolphram, "dados/c4_regra_188.data")
   call simulate_ACD(c4, 222, iterations_wolphram, L_wolphram, "dados/c4_regra_222.data")
   call simulate_ACD(c4, 250, iterations_wolphram, L_wolphram, "dados/c4_regra_250.data")


end program tarefa1




