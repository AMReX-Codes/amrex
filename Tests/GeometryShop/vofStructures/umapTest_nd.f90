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

  subroutine do_eb_work(lo, hi, bdmap, Nbd, ccmap, Ncc,&
    ktc, ktclo, ktchi, max_mv, ncomp) bind(C,name="do_eb_work")

    integer,     intent(in) :: lo(0:2), hi(0:2)
    integer,     intent(in) :: Nbd, Ncc
    type(ebbndrydata), intent(in) :: bdmap(0:Nbd-1)
    type(cutcell), intent(in) :: ccmap(0:Nbd-1)
    integer, intent(in) :: ktclo(0:2), ktchi(0:2), max_mv, ncomp
    integer, intent(in) :: ktc(ktclo(0):ktchi(0),ktclo(1):ktchi(1),ktclo(2):ktchi(2),0:max_mv-1,ncomp)

    integer :: i,j,k,n,L,key,d,tot,LL,RR,BB,TT
    real(amrex_real) :: vL, vR, vB, vT, vTOT

    do n = 1, ncomp
       do L = 0, max_mv-1
          do k = lo(2), hi(2)
             do j = lo(1), hi(1)
                do i = lo(0), hi(0)
                   key = ktc(i,j,k,L,n)
                   if (key .ge. 0) then

                      if (ccmap(key)%Nnbr(0,0).gt.0) then
                         LL = ccmap(key)%nbr(0,0,0)
                         if (LL .ge. 0) then
                            vL = bdmap(ktc(i-1,j,k,LL,n))%vol_frac
                         else
                            vL = 1.d0
                         endif
                      else
                         LL = -1
                         vL = 0.d0
                      endif

                      if (ccmap(key)%Nnbr(1,0).gt.0) then
                         RR = ccmap(key)%nbr(0,1,0)
                         if (RR .ge. 0) then
                            vR = bdmap(ktc(i+1,j,k,RR,n))%vol_frac
                         else
                            vR = 1.d0
                         endif
                      else
                         RR = -1
                         vR = 0.d0
                      endif

                      if (ccmap(key)%Nnbr(0,1).gt.0) then
                         BB = ccmap(key)%nbr(0,0,1)
                         if (BB .ge. 0) then
                            vB = bdmap(ktc(i,j-1,k,BB,n))%vol_frac
                         else
                            vB = 1.d0
                         endif
                      else
                         BB = -1
                         vB = 0.d0
                      endif

                      if (ccmap(key)%Nnbr(1,1).gt.0) then
                         TT = ccmap(key)%nbr(0,1,1)
                         if (TT .ge. 0) then
                            vT = bdmap(ktc(i,j+1,k,TT,n))%vol_frac
                         else
                            vT = 1.d0
                         endif
                      else
                         TT = -1
                         vT = 0.d0
                      endif

                      if (i.eq.37 .and. j.eq.57) then
                         print *,'VOL',vL,vR,vB,vT
                      endif

                   endif
                end do
             end do
          end do
       end do
    enddo


  end subroutine do_eb_work

end module umapTest_nd_module
