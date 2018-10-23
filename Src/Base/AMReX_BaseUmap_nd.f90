module baseumap_nd_module

use amrex_fort_module, only : amrex_real

  implicit none

contains


  function amrex_fort_umap_norm (lo, hi, src, src_sz, kt, ktlo, kthi, max_mv,ncomp, p) result(nrm) &
       bind(c,name='amrex_fort_umap_norm')
    integer, intent(in) :: lo(3), hi(3), ktlo(3), kthi(3), p, max_mv, ncomp, src_sz
    ! This dangerous - key_table_type is currently int, but could be changed to long
    ! mabye that should be an amrex_size_t?
    integer, intent(in) :: kt(ktlo(1):kthi(1),ktlo(2):kthi(2),ktlo(3):kthi(3),max_mv,ncomp)

    real(amrex_real), intent(in) :: src(0:src_sz-1) ! 0-indexing for c compatability
    real(amrex_real) :: nrm

    integer :: i,j,k,n, m, key

    nrm = 0.0_amrex_real
    if (.true.) then ! max norm - ignore p
       do n = 1, ncomp
           do m = 1, max_mv
              do       k = lo(3), hi(3)
                 do    j = lo(2), hi(2)
                    do i = lo(1), hi(1)
                       key = kt(i,j,k,m,n)
                       if (key .ge. 0) then
                           nrm = max(nrm, abs(src(key)))
                        endif
                    end do
                 end do
              end do
           end do
       enddo
    end if
  end function amrex_fort_umap_norm


  function amrex_fort_umap_norm_direct (lo, hi, src, src_sz, kt, ktlo, kthi, max_mv,ncomp, p) result(nrm) &
       bind(c,name='amrex_fort_umap_norm_direct')
    integer, intent(in) :: lo(3), hi(3), ktlo(3), kthi(3), p, max_mv, ncomp, src_sz
    ! This dangerous - key_table_type is currently int, but could be changed to long
    ! mabye that should be an amrex_size_t?
    integer, intent(in) :: kt(ktlo(1):kthi(1),ktlo(2):kthi(2),ktlo(3):kthi(3),max_mv,ncomp)

    real(amrex_real), intent(in) :: src(src_sz)
    real(amrex_real) :: nrm

    integer :: i,j,k,n, m, key

    nrm = 0.0_amrex_real
    if (.true.) then ! max norm - ignore p
        do m = 1, src_sz
           ! No way currently to check if this element is within the valid box
            nrm = max(nrm, abs(src(m)))
        enddo
    end if
  end function amrex_fort_umap_norm_direct
end module
