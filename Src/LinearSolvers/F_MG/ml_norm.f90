module ml_norm_module

  use multifab_module

  implicit none

  public :: ml_norm_inf, ml_norm_l2

contains

  function ml_norm_inf(mf, mask, local) result(r)
    type( multifab), intent(in)   :: mf(:)
    type(lmultifab), intent(in)   :: mask(:)
    logical, intent(in), optional :: local
    real(dp_t)                    :: r, r1
    integer                       :: n,nlevs
    nlevs = size(mf)
    r = norm_inf(mf(nlevs),local=local)
    do n = nlevs-1, 1, -1
       r1 = norm_inf(mf(n), mask(n), local=local)
       r = max(r1, r)
    end do
  end function ml_norm_inf

  function ml_norm_l2(mf, rr, mask) result(r)
    type( multifab), intent(in) :: mf(:)
    integer                     :: rr(:,:)
    type(lmultifab), intent(in) :: mask(:)
    real(dp_t)                  :: r
    integer                     :: n,nlevs
    nlevs = size(mf)
    r = norm_l2(mf(nlevs))**2
    do n = nlevs-1, 1, -1
       r =  r / product(rr(n,:)) &
          + norm_l2(mf(n), mask = mask(n))**2
    end do
    r = sqrt(r)
  end function ml_norm_l2

end module ml_norm_module
