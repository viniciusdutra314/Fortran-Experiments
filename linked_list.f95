module linked_list
   implicit none
   type :: particle_node
      integer :: position(2)
      type(particle_node), pointer :: next => null()
   end type particle_node

   type :: particle_list
      type(particle_node), pointer :: head => null()
      type(particle_node), pointer :: tail => null()
      integer :: size=0
   end type particle_list

contains
   subroutine add_particle_at_the_end(list, position)
      type(particle_list), intent(inout) :: list
      integer, intent(in) :: position(2)
      list%size = list%size + 1
      if (.not. associated(list%head)) then
         allocate(list%head)
         list%head%position = position
         list%head%next => null()
         list%tail => list%head
      else
         allocate(list%tail%next)
         list%tail%next%position = position
         list%tail%next%next => null()
         list%tail => list%tail%next
      end if
   end subroutine add_particle_at_the_end

   subroutine delete_particle(list, current_node, past_node)
      type(particle_list), intent(inout) :: list
      type(particle_node), pointer, intent(inout) :: current_node
      type(particle_node), pointer, intent(inout) :: past_node

      list%size = list%size - 1
      if (associated(current_node, list%head)) then
         list%head => current_node%next
         deallocate(current_node)
         current_node => list%head
      else
         past_node%next => current_node%next

         if (associated(current_node, list%tail)) then
            list%tail => past_node
         end if

         deallocate(current_node)
         current_node => past_node%next
      end if
   end subroutine delete_particle


   subroutine clear_list(list)
      type(particle_list), intent(inout) :: list
      type(particle_node), pointer :: current, next_node

      current => list%head
      do while (associated(current))
         next_node => current%next
         deallocate(current)
         current => next_node
      end do
      list%head => null()
      list%tail => null()
      list%size = 0
   end subroutine clear_list
end module linked_list

