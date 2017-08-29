subroutine amrex_compute_best_partition(cost, clo, chi, &
                                        lo, hi, total_cost, dir, &
                                        cost_left, cost_right, split) &
  bind(c,name='amrex_compute_best_partition')
    
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  integer              :: clo(3)
  integer              :: chi(3)    
  real(amrex_real)     :: cost(clo(1):chi(1), clo(2):chi(2), clo(3):chi(3))
  integer              :: lo(3)
  integer              :: hi(3)
  real(amrex_real), value :: total_cost
  integer, value          :: dir
  real(amrex_real)        :: cost_left
  real(amrex_real)        :: cost_right
  integer, intent(out)    :: split
  
  integer i,j,k
  real(amrex_real) target_cost, cost_sum
  target_cost = 0.5d0 * total_cost
  cost_sum = 0.d0
  
  if (dir .eq. 0) then
     
     split = lo(1)
     do i = lo(1), hi(1)-1
        do j = lo(2), hi(2)
           do k = lo(3), hi(3)
              cost_sum = cost_sum + cost(i, j, k)
           end do
        end do
        if (cost_sum .ge. target_cost) then
           split = i
           cost_left = cost_sum
           cost_right = total_cost - cost_left
           exit
        end if
     end do
     
  else if (dir .eq. 1) then
     
     split = lo(2)
     do j = lo(2), hi(2)-1
        do i = lo(1), hi(1)
           do k = lo(3), hi(3) 
              cost_sum = cost_sum + cost(i, j, k)
           end do
        end do
        if (cost_sum .ge. target_cost) then
           split = j
           cost_left = cost_sum
           cost_right = total_cost - cost_left
           exit
        end if
     end do
     
  else if (dir .eq. 2) then
     
     split = lo(3)
     do k = lo(3), hi(3)-1
        do j = lo(2), hi(2)
           do i = lo(1), hi(1) 
              cost_sum = cost_sum + cost(i, j, k)
           end do
        end do
        if (cost_sum .ge. target_cost) then
           split = k
           cost_left = cost_sum
           cost_right = total_cost - cost_left
           exit
        end if
     end do
     
  end if
  
end subroutine amrex_compute_best_partition
