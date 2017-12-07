subroutine amrex_compute_best_partition(cost, clo, chi, &
                                    lo, hi, total_cost, dir, &
                                    cost_left, cost_right, split) &
  bind(c,name='amrex_compute_best_partition')
  
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  integer              :: clo(1)
  integer              :: chi(1)    
  real(amrex_real)     :: cost(clo(1):chi(1))
  integer              :: lo(1)
  integer              :: hi(1)
  real(amrex_real), value :: total_cost
  integer, value          :: dir
  real(amrex_real)        :: cost_left
  real(amrex_real)        :: cost_right
  integer, intent(out)    :: split
  
  integer i
  real(amrex_real) target_cost, cost_sum
  target_cost = 0.5d0 * total_cost
  cost_sum = 0.d0
  
  split = lo(1)
  do i = lo(1), hi(1)
     cost_sum = cost_sum + cost(i)
     if (cost_sum .ge. target_cost) then
        split = i
        cost_left = cost_sum
        cost_right = total_cost - cost_left
        exit
     end if
  end do
  
end subroutine amrex_compute_best_partition

subroutine amrex_compute_cost(pcounts, cost, lo, hi, cell_weight) &
     bind(c,name='amrex_compute_cost')

  use iso_c_binding
  use amrex_fort_module, only : amrex_real

  integer                       :: lo(1)
  integer                       :: hi(1)
  real(amrex_real)              :: pcounts(lo(1):hi(1))
  real(amrex_real)              :: cost(lo(1):hi(1))
  real(amrex_real), value       :: cell_weight
  
  integer i
  
  do i = lo(1), hi(1)
     cost(i) = pcounts(i)*pcounts(i) + cell_weight
  end do
  
end subroutine amrex_compute_cost

subroutine amrex_set_box_cost(cost, clo, chi, &
                              lo, hi, box_cost) &

  bind(c,name='amrex_set_box_cost')
    
  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  integer              :: clo(1)
  integer              :: chi(1)    
  real(amrex_real)     :: cost(clo(1):chi(1))
  integer              :: lo(1)
  integer              :: hi(1)
  real(amrex_real)     :: box_cost
  
  integer i
  box_cost = 0.d0

  do i = lo(1), hi(1) 
     box_cost = box_cost + cost(i)
  end do

end subroutine amrex_set_box_cost
