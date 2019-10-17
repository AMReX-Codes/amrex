
module warpx_openbc_module

  implicit none

  integer, parameter :: idecomp = 0 ! 0=xyz, 1=xy, 2=yz, 3=xz, 4=x, 5=y, 6=z.
  integer, parameter :: igfflag = 1 ! =0 for ordinary 1/r Green function;
                                    ! =1 for integrated Green function

  integer, save :: gb_lo(3), gb_hi(3)
  integer, save :: lc_lo(3), lc_hi(3)

contains

    subroutine warpx_openbc_decompose(glo, ghi, lo, hi) bind(c,name='warpx_openbc_decompose')

    use mpi

    integer, intent(in) :: glo(3), ghi(3)
    integer, intent(out) :: lo(3), hi(3)
    integer :: myrank, mprocs, ierr, npx, npy, npz

    call MPI_Comm_size(MPI_COMM_WORLD, mprocs, ierr);
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    gb_lo = glo + 1  ! +1 because AMReX's domain index starts with 0
    gb_hi = ghi + 2  ! +2 to nodalize it

    call procgriddecomp(mprocs, gb_lo(1),gb_hi(1), &
         &                      gb_lo(2),gb_hi(2), &
         &                      gb_lo(3),gb_hi(3), &
         &                      idecomp, npx, npy, npz)

    call decompose(myrank,mprocs, gb_lo(1),gb_hi(1), &
         &                        gb_lo(2),gb_hi(2), &
         &                        gb_lo(3),gb_hi(3), &
         idecomp,                 lc_lo(1), lc_hi(1), &
         &                        lc_lo(2), lc_hi(2), &
         &                        lc_lo(3), lc_hi(3))

    lo = lc_lo - 1  ! AMReX's domain index starts with zero.
    hi = lc_hi - 1

  end subroutine warpx_openbc_decompose

  subroutine warpx_openbc_potential(rho, phi, dx) &
       bind(C, name="warpx_openbc_potential")

    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: rho(lc_lo(1):lc_hi(1),lc_lo(2):lc_hi(2),lc_lo(3):lc_hi(3))
    double precision, intent(out) :: phi(lc_lo(1):lc_hi(1),lc_lo(2):lc_hi(2),lc_lo(3):lc_hi(3))

    integer :: ierr

    call openbcpotential(rho, phi, dx(1), dx(2), dx(3), &
         lc_lo(1), lc_hi(1), lc_lo(2), lc_hi(2), lc_lo(3), lc_hi(3), &
         gb_lo(1), gb_hi(1), gb_lo(2), gb_hi(2), gb_lo(3), gb_hi(3), &
         idecomp, igfflag, ierr)

  end subroutine warpx_openbc_potential

end module warpx_openbc_module
