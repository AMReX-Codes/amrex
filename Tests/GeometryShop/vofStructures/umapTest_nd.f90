module umapTest_nd_module

  use amrex_fort_module, only : amrex_real, dim=>bl_spacedim
  use amrex_ebstruct_module, only : cutface, cutcell
  use umap_test_module, only : facedata, ebbndrydata

  implicit none

contains

  function fort_umap_norm (lo, hi, src, src_sz, kt, ktlo, kthi, max_mv,ncomp, p) result(nrm) &
       bind(c,name='fort_umap_norm')
    integer, intent(in) :: lo(3), hi(3), ktlo(3), kthi(3), p, max_mv, ncomp, src_sz
    integer, intent(in) :: kt(ktlo(1):kthi(1),ktlo(2):kthi(2),ktlo(3):kthi(3),max_mv,ncomp)

    real(amrex_real), intent(in) :: src(0:src_sz-1) ! 0-indexing for c compatability
    real(amrex_real) :: nrm

    integer :: i,j,k,n, m, key

    nrm = 0.0_amrex_real
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
  end function fort_umap_norm

  subroutine do_eb_work(lo, hi, cmap, ncmap, fmapI, fmapJ, nfmap, &
    ktc, ktclo, ktchi, max_mv, ncomp) bind(C,name="do_eb_work")

    integer,     intent(in) :: lo(0:2), hi(0:2)
    integer,     intent(in) :: ncmap, nfmap(0:2)
    type(ebbndrydata), intent(in) :: cmap(0:ncmap-1)
    type(facedata), intent(in) :: fmapI(0:nfmap(0)-1)
    type(facedata), intent(in) :: fmapJ(0:nfmap(1)-1)
    integer, intent(in) :: ktclo(0:2), ktchi(0:2), max_mv, ncomp
    integer, intent(in) :: ktc(ktclo(0):ktchi(0),ktclo(1):ktchi(1),ktclo(2):ktchi(2),max_mv,ncomp)

    integer :: i,j,k,n, m, key

    do n = 1, ncomp
       do m = 1, max_mv
          do k = lo(2), hi(2)
             do j = lo(1), hi(1)
                do i = lo(0), hi(0)
                   key = ktc(i,j,k,m,n)
                   if (key .ge. 0) then
                      print *,'found ',cmap(key)
                   endif
                end do
             end do
          end do
       end do
    enddo


  end subroutine do_eb_work

end module umapTest_nd_module
