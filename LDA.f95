module lda_module
   implicit none
contains
   subroutine initial_seed(aggregate)
      logical, target, intent(out) :: aggregate(:,:,:)
      integer :: L
      L=size(aggregate,1)
      aggregate = .false.
      if (size(aggregate,3) == 1) then
         aggregate(L/2, L/2, 1) = .true.
      else
         aggregate(L/2, L/2, L/2) = .true.
      end if
   end subroutine initial_seed

   function initial_position(R_inic, dimension,L) result(position)
      integer, intent(in) :: R_inic,dimension,L
      integer :: position(dimension)
      real :: phi,cos_theta,sin_theta
      real :: pi=3.14159265359
      if (dimension == 2) then
         call random_number(phi)
         phi=phi*2*pi
         position(1) = L/2+NINT(R_inic* cos(phi))
         position(2) = L/2+NINT(R_inic* sin(phi))
         position(3) = 1
      else if (dimension == 3) then
         call random_number(phi)
         call random_number(cos_theta)
         phi = phi * 2.0 * pi
         cos_theta = 2.0 * cos_theta - 1.0 ![-1,1]
         sin_theta=sqrt(1.0 - cos_theta**2)
         position(1) = L/2 + NINT(R_inic * cos(phi)*sin_theta)
         position(2) = L/2 + NINT(R_inic * sin(phi)*sin_theta)
         position(3) = L/2 + NINT(R_inic*cos_theta)
      end if
   end function initial_position

   subroutine random_step(particle,dimension)
      integer, intent(inout) :: particle(3)
      integer, intent(in) :: dimension
      real :: random
      integer :: random_sign,direction
      call random_number(random)
      if (random < 0.5) then
         random_sign = -1
      else
         random_sign = 1
      end if

      if (dimension == 2) then
         call random_number(random)
         if (random < 0.5) then
            direction =1
         else
            direction = 2
         end if
      else
         call random_number(random)
         if (random < 0.33333) then
            direction = 1
         else if (random < 0.66666) then
            direction = 2
         else
            direction = 3
         end if
      end if
      particle(direction) = particle(direction) + random_sign
   end subroutine random_step

   subroutine add_particle(aggregate, particle)
      logical, intent(inout) :: aggregate(:,:,:)
      integer, intent(in) :: particle(3)
      aggregate(particle(1),particle(2),particle(3)) = .true.
   end subroutine add_particle

   function get_distance(particle,L,dimension) result(distance)
      integer, intent(in) :: particle(3),L,dimension
      integer :: center(3)
      integer :: distance
      if (dimension == 2) then
         center=(/L/2,L/2,1/)
      else
         center=(/L/2,L/2,L/2/)
      end if
      distance=NINT(norm2(real(particle-center)))

   end function get_distance

   function particle_collided(particle, aggregate,d) result(collision)
      logical, intent(in) :: aggregate(:,:,:)
      integer, intent(in) :: particle(3),d
      integer, dimension(size(particle)) :: offset,possible_neighbor
      integer :: i,sign
      logical :: collision
      collision=.false.
      do i=1,d
         do sign=-1,1,2
            offset=0
            offset(i) = sign
            possible_neighbor=particle+offset
            if (aggregate(possible_neighbor(1),possible_neighbor(2),possible_neighbor(3))) then
               collision=.true.
               return
            end if
         end do
      end do
   end function particle_collided

   function particle_out_of_bounds(particle, aggregate,d) result(out_of_bounds)
      logical, intent(in) :: aggregate(:,:,:)
      integer, intent(in) :: particle(3),d
      integer :: i
      logical :: out_of_bounds
      out_of_bounds = .false.
      do i=1,d
         if (.not.((1<particle(i)) .and. particle(i)<size(aggregate,i))) then
            out_of_bounds = .true.
            return
         end if
      end do
   end function particle_out_of_bounds

   subroutine run_LDA(L,num_particles,d,filename)
      integer, intent(in) :: L,num_particles,d
      character(len=*), intent(in) :: filename
      integer :: n,distance,R_inic,max_distance,io_final_result
      logical, allocatable:: aggregate(:,:,:)
      integer :: particle(3)
      if (d==2) then
         allocate(aggregate(L,L,1))
      else if (d==3) then
         allocate(aggregate(L,L,L))
      end if
      max_distance=0
      R_inic=5+max_distance

      call initial_seed(aggregate)
      do n=0,num_particles-1
         particle=initial_position(R_inic,d,L)
         distance=get_distance(particle,L,d)
         particle_moving: do while(distance<=1.5*R_inic)
            call random_step(particle,d)
            distance=get_distance(particle,L,d)
            if (particle_out_of_bounds(particle,aggregate,d)) then
               exit particle_moving
            else if (particle_collided(particle,aggregate,d)) then
               call add_particle(aggregate,particle)
               if (distance>max_distance) then
                  max_distance=distance
                  R_inic=max_distance+20
               end if
               exit particle_moving
            end if
         end do particle_moving
      end do
      open(newunit=io_final_result, file=filename, status='replace')
      write(io_final_result,*) aggregate
      close(io_final_result)
   end subroutine run_LDA

end module lda_module



program LDA
   use lda_module
   implicit none
   integer :: L1,L2,num_particles1,num_particles2
   L1=1000
   L2=300

   num_particles1=100000
   num_particles2=num_particles1
   call run_LDA(L1,num_particles1,2,"aggregate2d.data")
   ! call run_LDA(L2,num_particles2,3,"aggregate3d.data")

end program LDA
