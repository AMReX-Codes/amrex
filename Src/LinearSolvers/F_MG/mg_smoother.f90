module mg_smoother_module

  use stencil_module
  use stencil_nodal_module
  use bl_types

  implicit none

  real(dp_t), private, parameter :: ZERO = 0.0_dp_t

  private dgtsl

contains

  subroutine gs_rb_smoother_1d(omega, ss, uu, ff, mm, lo, ng, n, skwd)
    integer, intent(in) :: lo(:)
    integer, intent(in) :: ng
    integer, intent(in) :: n
    real (kind = dp_t), intent(in)    :: omega
    real (kind = dp_t), intent(in)    :: ff(lo(1):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:)
    real (kind = dp_t), intent(in)    :: ss(lo(1):,0:)
    integer            ,intent(in)    :: mm(lo(1):)
    logical, intent(in), optional :: skwd

    real (kind = dp_t) dd
    integer ioff, i, hi(size(lo))
    integer, parameter ::  XBC = 3
    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    hi = ubound(uu)-ng

    if ( lskwd ) then
       ioff = 0; if ( mod(lo(1), 2) /= n ) ioff = 1
       do i = lo(1) + ioff, hi(1), 2
          dd = ss(i,0) * uu(i) + ss(i,1) * uu(i+1) + ss(i,2) * uu(i-1)
          if ( hi(1) > lo(1) ) then
             if (bc_skewed(mm(i),1,+1)) then
                dd = dd + ss(i,XBC)*uu(i+2)
             else if (bc_skewed(mm(i),1,-1)) then
                dd = dd + ss(i,XBC)*uu(i-2)
             end if
          end if
          if (abs(ss(i,0)) .gt. 0.0_dp_t) &
               uu(i) = uu(i) + omega/ss(i,0)*(ff(i) - dd)
       end do
    else
       ioff = 0; if ( mod(lo(1), 2) /= n ) ioff = 1
       do i = lo(1) + ioff, hi(1), 2
          dd = ss(i,0) * uu(i) + ss(i,1) * uu(i+1) + ss(i,2) * uu(i-1)
          if (abs(ss(i,0)) .gt. 0.0_dp_t) &
               uu(i) = uu(i) + omega/ss(i,0)*(ff(i) - dd)
       end do
    end if
  end subroutine gs_rb_smoother_1d

  subroutine gs_rb_smoother_2d(omega, ss, uu, ff, mm, lo, ng, n, skwd)
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    integer, intent(in) :: n
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(lo(1):, lo(2):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:)
    real (kind = dp_t), intent(in) :: ss(lo(1):, lo(2):, 0:)
    integer            ,intent(in) :: mm(lo(1):,lo(2):)
    logical, intent(in), optional :: skwd
    integer j, i
    integer hi(size(lo))
    integer ioff
    integer, parameter ::  XBC = 5, YBC = 6
    real (kind = dp_t) dd

    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    hi = ubound(uu)-ng

    if ( lskwd ) then
    !$OMP PARALLEL DO PRIVATE(j,i,ioff,dd)
       do j = lo(2),hi(2)
          ioff = 0; if ( mod(lo(1) + j, 2) /= n ) ioff = 1
          do i = lo(1) + ioff, hi(1), 2

             dd = ss(i,j,0)*uu(i,j) &
                  + ss(i,j,1) * uu(i+1,j) + ss(i,j,2) * uu(i-1,j) &
                  + ss(i,j,3) * uu(i,j+1) + ss(i,j,4) * uu(i,j-1)
             if (hi(1) > lo(1)) then
                if (bc_skewed(mm(i,j),1,+1)) then
                   dd = dd + ss(i,j,XBC)*uu(i+2,j)
                else if (bc_skewed(mm(i,j),1,-1)) then
                   dd = dd + ss(i,j,XBC)*uu(i-2,j)
                end if
             end if
             if (hi(2) > lo(2)) then
                if (bc_skewed(mm(i,j),2,+1)) then
                   dd = dd + ss(i,j,YBC)*uu(i,j+2)
                else if (bc_skewed(mm(i,j),2,-1)) then
                   dd = dd + ss(i,j,YBC)*uu(i,j-2)
                end if
             end if

             if (abs(ss(i,j,0)) .gt. 0.0_dp_t) &
                  uu(i,j) = uu(i,j) + omega/ss(i,j,0)*(ff(i,j) - dd)
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(j,i,ioff,dd)
       do j = lo(2),hi(2)
          ioff = 0; if ( mod(lo(1) + j, 2) /= n ) ioff = 1
          do i = lo(1) + ioff, hi(1), 2

             dd = ss(i,j,0)*uu(i,j) &
                  + ss(i,j,1) * uu(i+1,j) + ss(i,j,2) * uu(i-1,j) &
                  + ss(i,j,3) * uu(i,j+1) + ss(i,j,4) * uu(i,j-1)
             if (abs(ss(i,j,0)) .gt. 0.0_dp_t) &
                  uu(i,j) = uu(i,j) + omega/ss(i,j,0)*(ff(i,j) - dd)
          end do
       end do
    !$OMP END PARALLEL DO
    end if

  end subroutine gs_rb_smoother_2d

  subroutine gs_rb_smoother_3d(omega, ss, uu, ff, mm, lo, ng, n, skwd)
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    integer, intent(in) :: n
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(lo(1):,lo(2):,lo(3):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in) :: ss(lo(1):, lo(2):, lo(3):, 0:)
    integer            ,intent(in) :: mm(lo(1):,lo(2):,lo(3):)
    logical, intent(in), optional :: skwd
    integer i, j, k, ioff
    integer hi(size(lo))
    real (kind = dp_t) dd
    integer, parameter ::  XBC = 7, YBC = 8, ZBC = 9
    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    hi = ubound(uu)-ng

    if ( lskwd ) then
       !$OMP PARALLEL DO PRIVATE(k,j,i,ioff,dd)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             ioff = 0; if ( mod (lo(1) + j + k, 2) /= n ) ioff = 1
             do i = lo(1)+ioff, hi(1), 2

                dd = ss(i,j,k,0)*uu(i,j,k)
                dd = dd + ss(i,j,k,1)*uu(i+1,j,k) + ss(i,j,k,2)*uu(i-1,j,k)
                dd = dd + ss(i,j,k,3)*uu(i,j+1,k) + ss(i,j,k,4)*uu(i,j-1,k)
                dd = dd + ss(i,j,k,5)*uu(i,j,k+1) + ss(i,j,k,6)*uu(i,j,k-1)
                if (hi(1) > lo(1)) then
                   if (bc_skewed(mm(i,j,k),1,+1)) then
                      dd = dd + ss(i,j,k,XBC)*uu(i+2,j,k)
                   else if (bc_skewed(mm(i,j,k),1,-1)) then
                      dd = dd + ss(i,j,k,XBC)*uu(i-2,j,k)
                   end if
                end if

                if (hi(2) > lo(2)) then
                   if (bc_skewed(mm(i,j,k),2,+1)) then
                      dd = dd + ss(i,j,k,YBC)*uu(i,j+2,k)
                   else if (bc_skewed(mm(i,j,k),2,-1)) then
                      dd = dd + ss(i,j,k,YBC)*uu(i,j-2,k)
                   end if
                end if

                if (hi(3) > lo(3)) then
                   if (bc_skewed(mm(i,j,k),3,+1)) then
                      dd = dd + ss(i,j,k,ZBC)*uu(i,j,k+2)
                   else if (bc_skewed(mm(i,j,k),3,-1)) then
                      dd = dd + ss(i,j,k,ZBC)*uu(i,j,k-2)
                   end if
                end if
                if (abs(ss(i,j,k,0)) .gt. 0.0_dp_t) &
                     uu(i,j,k) = uu(i,j,k) + omega/ss(i,j,k,0)*(ff(i,j,k) - dd)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(k,j,i,ioff,dd)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             ioff = 0; if ( mod (lo(1) + j + k, 2) /= n ) ioff = 1
             do i = lo(1)+ioff, hi(1), 2

                dd = ss(i,j,k,0)*uu(i,j,k)
                dd = dd + ss(i,j,k,1)*uu(i+1,j,k) + ss(i,j,k,2)*uu(i-1,j,k)
                dd = dd + ss(i,j,k,3)*uu(i,j+1,k) + ss(i,j,k,4)*uu(i,j-1,k)
                dd = dd + ss(i,j,k,5)*uu(i,j,k+1) + ss(i,j,k,6)*uu(i,j,k-1)
                if (abs(ss(i,j,k,0)) .gt. 0.0_dp_t) &
                     uu(i,j,k) = uu(i,j,k) + omega/ss(i,j,k,0)*(ff(i,j,k) - dd)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine gs_rb_smoother_3d

  subroutine nodal_smoother_1d(omega, ss, uu, ff, mm, lo, ng)
    integer, intent(in) :: lo(:)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)    :: omega
    real (kind = dp_t), intent(in)    :: ff(lo(1)-1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:)
    real (kind = dp_t), intent(in)    :: ss(lo(1):,0:)
    integer            ,intent(in)    :: mm(lo(1):)

    real (kind = dp_t) dd
    integer i, hi(size(lo))

    hi = ubound(uu)-ng

    !   Leave Dirichlet boundaries untouched.
    i = lo(1)
    if (.not. bc_dirichlet(mm(i),1,0)) then
       dd = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
       uu(i) = uu(i) + omega/ss(i,0)*(ff(i) - dd)
    end if

    do i = lo(1)+1,hi(1)-1
       dd = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
       uu(i) = uu(i) + omega/ss(i,0)*(ff(i) - dd)
    end do

    i = hi(1)
    !   Leave Dirichlet boundaries untouched.
    if (.not. bc_dirichlet(mm(i),1,0)) then
       dd = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
       uu(i) = uu(i) + omega/ss(i,0)*(ff(i) - dd)
    end if

  end subroutine nodal_smoother_1d

  subroutine nodal_smoother_2d(omega, ss, uu, ff, mm, lo, ng)
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(lo(1)-1:, lo(2)-1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:)
    real (kind = dp_t), intent(in) :: ss(lo(1):, lo(2):, 0:)
    integer            ,intent(in) :: mm(lo(1):,lo(2):)
    integer j, i
    integer hi(size(lo))
    real (kind = dp_t) dd

    hi = ubound(uu)-ng
    dd = ZERO

    call impose_neumann_bcs_2d(uu,mm,lo,ng)

    !$OMP PARALLEL DO PRIVATE(j,i,dd)
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          if (.not. bc_dirichlet(mm(i,j),1,0)) then
             dd = ss(i,j,0)*uu(i,j) + ss(i,j,1) * uu(i-1,j-1) &
                  + ss(i,j,2) * uu(i  ,j-1) &
                  + ss(i,j,3) * uu(i+1,j-1) &
                  + ss(i,j,4) * uu(i-1,j  ) &
                  + ss(i,j,5) * uu(i+1,j  ) &
                  + ss(i,j,6) * uu(i-1,j+1) &
                  + ss(i,j,7) * uu(i  ,j+1) &
                  + ss(i,j,8) * uu(i+1,j+1)
             uu(i,j) = uu(i,j) + omega/ss(i,j,0)*(ff(i,j) - dd)
             !         write(6,1000) i,j,uu(i,j),ff(i,j),1./ss(i,j,0)
          end if
       end do
    end do
    !$OMP END PARALLEL DO
    ! 1000 format('NEW U ',i2,1x,i2,1x,f14.9,1x,f14.9,1x,f14.9)

  end subroutine nodal_smoother_2d

  subroutine nodal_smoother_3d(omega, ss, uu, ff, mm, lo, ng)
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in) :: ss(lo(1):, lo(2):, lo(3):, 0:)
    integer            ,intent(in) :: mm(lo(1):,lo(2):,lo(3):)
    integer i, j, k
    integer hi(size(lo))
    real (kind = dp_t) dd

    hi(1) = lo(1) + size(mm,dim=1)-1
    hi(2) = lo(2) + size(mm,dim=2)-1
    hi(3) = lo(3) + size(mm,dim=3)-1

    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
                dd = ss(i,j,k,0)*uu(i,j,k) &
                     + ss(i,j,k, 1) * uu(i-1,j-1,k-1) + ss(i,j,k, 2) * uu(i  ,j-1,k-1) &
                     + ss(i,j,k, 3) * uu(i+1,j-1,k-1) + ss(i,j,k, 4) * uu(i-1,j  ,k-1) &
                     + ss(i,j,k, 5) * uu(i+1,j  ,k-1) + ss(i,j,k, 6) * uu(i-1,j+1,k-1) &
                     + ss(i,j,k, 7) * uu(i  ,j+1,k-1) + ss(i,j,k, 8) * uu(i+1,j+1,k-1) &
                     + ss(i,j,k, 9) * uu(i-1,j-1,k  ) + ss(i,j,k,10) * uu(i+1,j-1,k  ) &
                     + ss(i,j,k,11) * uu(i-1,j+1,k  ) + ss(i,j,k,12) * uu(i+1,j+1,k  ) &
                     + ss(i,j,k,13) * uu(i-1,j-1,k+1) + ss(i,j,k,14) * uu(i  ,j-1,k+1) &
                     + ss(i,j,k,15) * uu(i+1,j-1,k+1) + ss(i,j,k,16) * uu(i-1,j  ,k+1) &
                     + ss(i,j,k,17) * uu(i+1,j  ,k+1) + ss(i,j,k,18) * uu(i-1,j+1,k+1) &
                     + ss(i,j,k,19) * uu(i  ,j+1,k+1) + ss(i,j,k,20) * uu(i+1,j+1,k+1) 

                !         Add faces (only non-zero for non-uniform dx)
                dd = dd + &
                     ss(i,j,k,21) * uu(i-1,j  ,k  ) + ss(i,j,k,22) * uu(i+1,j  ,k  ) &
                     + ss(i,j,k,23) * uu(i  ,j-1,k  ) + ss(i,j,k,24) * uu(i  ,j+1,k  ) &
                     + ss(i,j,k,25) * uu(i  ,j  ,k-1) + ss(i,j,k,26) * uu(i  ,j  ,k+1)

                uu(i,j,k) = uu(i,j,k) + omega/ss(i,j,k,0)*(ff(i,j,k) - dd)
                !         write(6,1000) i,j,k,uu(i,j,k),ff(i,j,k),1./ss(i,j,k,0)
             end if
          end do
       end do
       !    print *,' '
    end do

    ! 1000 format('NEW U ',i2,1x,i2,1x,i2,1x,f14.9,1x,f14.9,1x,f14.9)

  end subroutine nodal_smoother_3d

  subroutine gs_lex_smoother_1d(omega, ss, uu, ff, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)    :: omega
    real (kind = dp_t), intent(in)    :: ff(:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:)
    real (kind = dp_t), intent(in)    :: ss(:,0:)
    integer            ,intent(in)    :: mm(:)
    logical, intent(in), optional :: skwd
    real (kind = dp_t) dd
    integer i
    logical :: lskwd
    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    do i = 1, size(ff,dim=1)
       dd = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
       uu(i) = uu(i) + omega/ss(i,0)*(ff(i) - dd)
    end do

  end subroutine gs_lex_smoother_1d

  subroutine gs_lex_smoother_2d(omega, ss, uu, ff, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(:, :)
    real (kind = dp_t), intent(inout) :: uu(1-ng:, 1-ng:)
    real (kind = dp_t), intent(in) :: ss(:, :, 0:)
    integer            ,intent(in) :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer j, i
    real (kind = dp_t) dd
    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    !$OMP PARALLEL DO PRIVATE(j,i,dd)
    do j = 1, size(ff,dim=2)
       do i = 1, size(ff, dim=1)
          dd = &
               + ss(i,j,0)*uu(i  ,j  ) &
               + ss(i,j,1)*uu(i+1,j  ) &
               + ss(i,j,2)*uu(i-1,j  ) &
               + ss(i,j,3)*uu(i  ,j+1) &
               + ss(i,j,4)*uu(i  ,j-1)
          uu(i,j) = uu(i,j) + omega/ss(i,j,0)*(ff(i,j) - dd)
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine gs_lex_smoother_2d

  subroutine gs_lex_smoother_3d(omega, ss, uu, ff, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(:,:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(in) :: ss(:,:,:,0:)
    integer            ,intent(in) :: mm(:,:,:)
    logical, intent(in), optional :: skwd
    integer i, j, k
    real (kind = dp_t) dd
    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    !$OMP PARALLEL DO PRIVATE(k,j,i,dd)
    do k = 1, size(ff,dim=3)
       do j = 1, size(ff,dim=2)
          do i = 1, size(ff,dim=1)
                dd = &
                     + ss(i,j,k, 0)*uu(  i,j  ,k  ) &
                     + ss(i,j,k, 1)*uu(i+1,j  ,k  ) &
                     + ss(i,j,k, 2)*uu(i-1,j  ,k  ) &
                     + ss(i,j,k, 3)*uu(i  ,j+1,k  ) &
                     + ss(i,j,k, 4)*uu(i  ,j-1,k  ) &
                     + ss(i,j,k, 5)*uu(i  ,j  ,k-1) &
                     + ss(i,j,k, 6)*uu(i  ,j  ,k+1) 
                uu(i,j,k) = uu(i,j,k) + omega/ss(i,j,k,0)*(ff(i,j,k) - dd)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine gs_lex_smoother_3d  

subroutine gs_lex_dense_smoother_1d(omega, ss, uu, ff, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)    :: omega
    real (kind = dp_t), intent(in)    :: ff(:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:)
    real (kind = dp_t), intent(in)    :: ss(:,0:)
    integer            ,intent(in)    :: mm(:)
    logical, intent(in), optional :: skwd
    real (kind = dp_t) dd
    integer i
    logical :: lskwd
    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    do i = 1, size(ff,dim=1)
       dd = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
       uu(i) = uu(i) + omega/ss(i,0)*(ff(i) - dd)
    end do


  end subroutine gs_lex_dense_smoother_1d

  subroutine gs_lex_dense_smoother_2d(omega, ss, uu, ff, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(:, :)
    real (kind = dp_t), intent(inout) :: uu(1-ng:, 1-ng:)
    real (kind = dp_t), intent(in) :: ss(:, :, 0:)
    integer            ,intent(in) :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer j, i
    real (kind = dp_t) dd
    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    !$OMP PARALLEL DO PRIVATE(j,i,dd)
    do j = 1, size(ff,dim=2)
       do i = 1, size(ff, dim=1)
          dd = &
               + ss(i,j,1)*uu(i-1,j-1) &
               + ss(i,j,2)*uu(i  ,j-1) &
               + ss(i,j,3)*uu(i+1,j-1) &

               + ss(i,j,4)*uu(i-1,j  ) &
               + ss(i,j,0)*uu(i  ,j  ) &
               + ss(i,j,5)*uu(i+1,j  ) &

               + ss(i,j,6)*uu(i-1,j+1) &
               + ss(i,j,7)*uu(i  ,j+1) &
               + ss(i,j,8)*uu(i+1,j+1)
          uu(i,j) = uu(i,j) + omega/ss(i,j,0)*(ff(i,j) - dd)
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine gs_lex_dense_smoother_2d

  subroutine gs_lex_dense_smoother_3d(omega, ss, uu, ff, mm, ng, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ff(:,:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(in) :: ss(:,:,:,0:)
    integer            ,intent(in) :: mm(:,:,:)
    logical, intent(in), optional :: skwd
    integer i, j, k
    real (kind = dp_t) dd
    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    !$OMP PARALLEL DO PRIVATE(k,j,i,dd)
    do k = 1, size(ff,dim=3)
       do j = 1, size(ff,dim=2)
          do i = 1, size(ff,dim=1)
                dd = &
                     + ss(i,j,k, 0)*uu(i  ,j  ,k  ) &
                     + ss(i,j,k, 1)*uu(i-1,j-1,k-1) &
                     + ss(i,j,k, 2)*uu(i  ,j-1,k-1) &
                     + ss(i,j,k, 3)*uu(i+1,j-1,k-1) &
                     + ss(i,j,k, 4)*uu(i-1,j  ,k-1) &
                     + ss(i,j,k, 5)*uu(i  ,j  ,k-1) &
                     + ss(i,j,k, 6)*uu(i+1,j  ,k-1) &
                     + ss(i,j,k, 7)*uu(i-1,j+1,k-1) &
                     + ss(i,j,k, 8)*uu(i  ,j+1,k-1) &
                     + ss(i,j,k, 9)*uu(i+1,j+1,k-1) &
                     
                     + ss(i,j,k,10)*uu(i-1,j-1,k  ) &
                     + ss(i,j,k,11)*uu(i  ,j-1,k  ) &
                     + ss(i,j,k,12)*uu(i+1,j-1,k  ) &
                     + ss(i,j,k,13)*uu(i-1,j  ,k  ) &
                     + ss(i,j,k, 0)*uu(i  ,j  ,k  ) &
                     + ss(i,j,k,14)*uu(i+1,j  ,k  ) &
                     + ss(i,j,k,15)*uu(i  ,j+1,k  ) &
                     + ss(i,j,k,16)*uu(i-1,j+1,k  ) &
                     + ss(i,j,k,17)*uu(i+1,j+1,k  ) &
         
                     + ss(i,j,k,18)*uu(i-1,j-1,k+1) &
                     + ss(i,j,k,19)*uu(i  ,j-1,k+1) &
                     + ss(i,j,k,20)*uu(i+1,j-1,k+1) &
                     + ss(i,j,k,21)*uu(i-1,j  ,k+1) &
                     + ss(i,j,k,22)*uu(i  ,j  ,k+1) &
                     + ss(i,j,k,23)*uu(i+1,j  ,k+1) &
                     + ss(i,j,k,24)*uu(i-1,j+1,k+1) &
                     + ss(i,j,k,25)*uu(i  ,j+1,k+1) &
                     + ss(i,j,k,26)*uu(i+1,j+1,k+1) 

                uu(i,j,k) = uu(i,j,k) + omega/ss(i,j,k,0)*(ff(i,j,k) - dd)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine gs_lex_dense_smoother_3d

  subroutine jac_smoother_1d(omega, ss, uu, ff, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) ::  omega
    real (kind = dp_t), intent(in) ::  ss(:,0:)
    real (kind = dp_t), intent(in) ::  ff(:)
    integer           , intent(in) ::  mm(:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:)
    real (kind = dp_t):: wrk(size(ff,1))
    integer nx
    integer i
    real (kind = dp_t) :: dd
    integer, parameter ::  XBC = 3

    nx = size(ff,dim=1)
    do i = 1, nx
       dd = + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1) 
       if ( nx > 1 ) then
          if ( bc_skewed(mm(i),1,+1) ) then
             dd = dd + ss(i,XBC)*uu(i+2)
          else if ( bc_skewed(mm(i),1,-1) ) then
             dd = dd + ss(i,XBC)*uu(i-2)
          end if
       end if
       wrk(i) = ff(i) - dd
    end do
    do i = 1, nx
       uu(i) = uu(i) + omega*(wrk(i)/ss(i,0)-uu(i))
    end do

  end subroutine jac_smoother_1d

  subroutine jac_dense_smoother_1d(omega, ss, uu, ff, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) ::  omega
    real (kind = dp_t), intent(in) ::  ss(:,:)
    real (kind = dp_t), intent(in) ::  ff(:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:)
    real (kind = dp_t):: wrk(size(ff,1))
    integer nx, i

    nx = size(ff,dim=1)
    do i = 1, nx
       wrk(i) = ( ff(i) &
            - ss(i,1)*uu(i+1) - ss(i,3)*uu(i-1) &
            )
    end do
    do i = 1, nx
       uu(i) = uu(i) + omega*(wrk(i)/ss(i,2)-uu(i))
    end do

  end subroutine jac_dense_smoother_1d

  subroutine jac_smoother_2d(omega, ss, uu, ff, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) ::  omega
    real (kind = dp_t), intent(in) ::  ss(:,:,0:)
    real (kind = dp_t), intent(in) ::  ff(:,:)
    integer           , intent(in) ::  mm(:,:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:,1-ng:)
    real (kind = dp_t) :: wrk(size(ff,1),size(ff,2))
    integer nx, ny
    integer i, j
    real (kind = dp_t) :: dd
    integer, parameter ::  XBC = 5, YBC = 6

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)

    !$OMP PARALLEL PRIVATE(j,i,dd)
    !$OMP DO
    do j = 1, ny
       do i = 1, nx
          dd =    + ss(i,j,1) * uu(i+1,j) + ss(i,j,2) * uu(i-1,j)
          dd = dd + ss(i,j,3) * uu(i,j+1) + ss(i,j,4) * uu(i,j-1)
          if ( nx > 1 ) then
             if (bc_skewed(mm(i,j),1,+1)) then
                dd = dd + ss(i,j,XBC)*uu(i+2,j)
             else if (bc_skewed(mm(i,j),1,-1)) then
                dd = dd + ss(i,j,XBC)*uu(i-2,j)
             end if
          end if
          if ( ny > 1 ) then
             if (bc_skewed(mm(i,j),2,+1)) then
                dd = dd + ss(i,j,YBC)*uu(i,j+2)
             else if (bc_skewed(mm(i,j),2,-1)) then
                dd = dd + ss(i,j,YBC)*uu(i,j-2)
             end if
          end if
          wrk(i,j) = ff(i,j) - dd
       end do
    end do
    !$OMP END DO
    !$OMP DO
    do j = 1, ny
       do i = 1, nx
          uu(i,j) = uu(i,j) + omega*(wrk(i,j)/ss(i,j,0)-uu(i,j))
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine jac_smoother_2d

  subroutine jac_dense_smoother_2d(omega, ss, uu, ff, ng)
    integer nx, ny
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ss(:,:,0:)
    real (kind = dp_t), intent(in) :: ff(:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:)
    real (kind = dp_t) :: wrk(size(ff,1),size(ff,2))
    integer i, j

    nx=size(ss,dim=1)
    ny=size(ss,dim=2)
    !$OMP PARALLEL PRIVATE(j,i)
    !$OMP DO
    do j = 1, ny
       do i = 1, nx
          wrk(i,j) = ff(i,j) - (&
               + ss(i,j,1)*uu(i-1,j-1) &
               + ss(i,j,2)*uu(i  ,j-1) &
               + ss(i,j,3)*uu(i+1,j-1) &

               + ss(i,j,4)*uu(i-1,j  ) &
               + ss(i,j,5)*uu(i+1,j  ) &

               + ss(i,j,6)*uu(i-1,j+1) &
               + ss(i,j,7)*uu(i  ,j+1) &
               + ss(i,j,8)*uu(i+1,j+1) &
               )
       end do
    end do
    !$OMP END DO
    !$OMP DO
    do j = 1, ny
       do i = 1, nx
          uu(i,j) = uu(i,j) + omega*(wrk(i,j)/ss(i,j,0)-uu(i,j))
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine jac_dense_smoother_2d

  subroutine jac_smoother_3d(omega, ss, uu, ff, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(in) :: ff(:,:,:)
    integer           , intent(in) :: mm(:,:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), allocatable :: wrk(:,:,:)
    real (kind = dp_t) :: dd
    integer nx, ny, nz
    integer i, j, k
    integer, parameter ::  XBC = 7, YBC = 8, ZBC = 9

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)
    nz = size(ff,dim=3)
    allocate(wrk(nx,ny,nz))
    !$OMP PARALLEL PRIVATE(j,i,k,dd)
    !$OMP DO
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             dd =    + ss(i,j,k,1)*uu(i+1,j,k) + ss(i,j,k,2)*uu(i-1,j,k)
             dd = dd + ss(i,j,k,3)*uu(i,j+1,k) + ss(i,j,k,4)*uu(i,j-1,k)
             dd = dd + ss(i,j,k,5)*uu(i,j,k+1) + ss(i,j,k,6)*uu(i,j,k-1)

             if ( nx > 1 ) then
                if (bc_skewed(mm(i,j,k),1,+1)) then
                   dd = dd + ss(i,j,k,XBC)*uu(i+2,j,k)
                else if (bc_skewed(mm(i,j,k),1,-1)) then
                   dd = dd + ss(i,j,k,XBC)*uu(i-2,j,k)
                end if
             end if

             if ( ny > 1 ) then
                if (bc_skewed(mm(i,j,k),2,+1)) then
                   dd = dd + ss(i,j,k,YBC)*uu(i,j+2,k)
                else if (bc_skewed(mm(i,j,k),2,-1)) then
                   dd = dd + ss(i,j,k,YBC)*uu(i,j-2,k)
                end if
             end if

             if ( nz > 1 ) then
                if (bc_skewed(mm(i,j,k),3,+1)) then
                   dd = dd + ss(i,j,k,ZBC)*uu(i,j,k+2)
                else if (bc_skewed(mm(i,j,k),3,-1)) then
                   dd = dd + ss(i,j,k,ZBC)*uu(i,j,k-2)
                end if
             end if
             wrk(i,j,k) = ff(i,j,k) - dd
          end do
       end do
    end do
    !$OMP END DO
    !$OMP DO
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             uu(i,j,k) = uu(i,j,k) + omega*(wrk(i,j,k)/ss(i,j,k,0) - uu(i,j,k))
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine jac_smoother_3d

  subroutine jac_dense_smoother_3d(omega, ss, uu, ff, ng)
    integer nx, ny, nz
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(in) :: ff(:,:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), allocatable :: wrk(:,:,:)
    integer i, j, k

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)
    nz = size(ff,dim=3)
    allocate(wrk(nx,ny,nz))
    !$OMP PARALLEL PRIVATE(j,i,k)
    !$OMP DO
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             wrk(i,j,k) = ff(i,j,k) - ( &
                  + ss(i,j,k, 0)*uu(i  ,j  ,k)  &
                  + ss(i,j,k, 1)*uu(i-1,j-1,k-1) &
                  + ss(i,j,k, 2)*uu(i  ,j-1,k-1) &
                  + ss(i,j,k, 3)*uu(i+1,j-1,k-1) &
                  + ss(i,j,k, 4)*uu(i-1,j  ,k-1) &
                  + ss(i,j,k, 5)*uu(i  ,j  ,k-1) &
                  + ss(i,j,k, 6)*uu(i+1,j  ,k-1) &
                  + ss(i,j,k, 7)*uu(i-1,j+1,k-1) &
                  + ss(i,j,k, 8)*uu(i  ,j+1,k-1) &
                  + ss(i,j,k, 9)*uu(i+1,j+1,k-1) &
                  
                  + ss(i,j,k,10)*uu(i-1,j-1,k  ) &
                  + ss(i,j,k,11)*uu(i  ,j-1,k  ) &
                  + ss(i,j,k,12)*uu(i+1,j-1,k  ) &
                  + ss(i,j,k,13)*uu(i-1,j  ,k  ) &
                  + ss(i,j,k,14)*uu(i+1,j  ,k  ) &
                  + ss(i,j,k,15)*uu(i-1,j+1,k  ) &
                  + ss(i,j,k,16)*uu(i  ,j+1,k  ) &
                  + ss(i,j,k,17)*uu(i+1,j+1,k  ) &
                  
                  + ss(i,j,k,18)*uu(i-1,j-1,k+1) &
                  + ss(i,j,k,19)*uu(i  ,j-1,k+1) &
                  + ss(i,j,k,20)*uu(i+1,j-1,k+1) &
                  + ss(i,j,k,21)*uu(i-1,j  ,k+1) &
                  + ss(i,j,k,22)*uu(i  ,j  ,k+1) &
                  + ss(i,j,k,23)*uu(i+1,j  ,k+1) &
                  + ss(i,j,k,24)*uu(i-1,j+1,k+1) &
                  + ss(i,j,k,25)*uu(i  ,j+1,k+1) &
                  + ss(i,j,k,26)*uu(i+1,j+1,k+1) &
                  )
          end do
       end do
    end do
    !$OMP END DO
    !$OMP DO
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             uu(i,j,k) = uu(i,j,k) + omega*(wrk(i,j,k)/ss(i,j,k,14)-uu(i,j,k))
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine jac_dense_smoother_3d

  subroutine lgs_x_2d(omega, ss, uu, ff, lo, ng)
    integer, intent(in) :: lo(:), ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ss(lo(1):,lo(2):,0:)
    real (kind = dp_t), intent(in) :: ff(lo(1):,lo(2):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:)

    real (kind = dp_t) :: wrk(size(ff,dim=1),4)
    integer i, j, hi(size(lo))

    hi = ubound(uu)-ng

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          wrk(i,1:3) = ss(i,j,(/0,2,1/))
          wrk(i,4) = ff(i,j) - (ss(i,j,3)*uu(i,j+1) + ss(i,j,4)*uu(i,j-1))
       end do
       call dgtsl(wrk(:,2), wrk(:,1), wrk(:,3), wrk(:,4))
       do i = lo(1), hi(1)
          uu(i,j) = uu(i,j) + omega*(wrk(i,4)-uu(i,j))
       end do
    end do
  end subroutine lgs_x_2d

  subroutine lgs_y_2d(omega, ss, uu, ff, lo, ng)
    integer, intent(in) :: lo(:), ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ss(lo(1):,lo(2):,0:)
    real (kind = dp_t), intent(in) :: ff(lo(1):,lo(2):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t) :: wrk(size(ff,dim=2),4)
    integer i, j, hi(size(lo))
    hi = ubound(uu)-ng
    do i = lo(1), hi(1)
       do j = lo(2), hi(2)
          wrk(j,1:3) = ss(i,j,(/0,4,3/))
          wrk(j,4) = ff(i,j) -(ss(i,j,1)*uu(i+1,j) + ss(i,j,2)*uu(i-1,j))
       end do
       call dgtsl(wrk(:,2), wrk(:,1), wrk(:,3), wrk(:,4))
       do j = lo(2), hi(2)
          uu(i,j) = uu(i,j) - omega*(wrk(j,4)-uu(i,j))
       end do
    end do
  end subroutine lgs_y_2d

  subroutine lgs_rb_x_2d(omega, ss, uu, ff, lo, ng)
    integer, intent(in) :: lo(:), ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ss(lo(1):,lo(2):,0:)
    real (kind = dp_t), intent(in) :: ff(lo(1):,lo(2):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t) :: wrk(size(ff,dim=1),4)
    integer i, j, hi(size(lo))

    hi = ubound(uu)-ng

    ! Note: it is wasteful to use the 4 component of wrk

    do j = lo(2), hi(2), 2
       do i = lo(1), hi(1)
          wrk(i,1:3) = ss(i,j,(/0,2,1/))
          wrk(i,4) = ff(i,j) -(ss(i,j,3)*uu(i,j+1) + ss(i,j,4)*uu(i,j-1))
       end do
       call dgtsl(wrk(:,2), wrk(:,1), wrk(:,3), wrk(:,4))
       do i = lo(1), hi(1)
          uu(i,j) = uu(i,j) + omega*(wrk(i,4)-uu(i,j))
       end do
    end do

    do j = lo(2)+1, hi(2), 2
       do i = lo(1), hi(1)
          wrk(i,1:3) = ss(i,j,(/0,2,1/))
          wrk(i,4) = -(ss(i,j,3)*uu(i,j+1) + ss(i,j,4)*uu(i,j-1))
       end do
       call dgtsl(wrk(:,2), wrk(:,1), wrk(:,3), wrk(:,4))
       do i = lo(1), hi(1)
          uu(i,j) = uu(i,j) + omega*(wrk(i,4)-uu(i,j))
       end do
    end do
  end subroutine lgs_rb_x_2d

  subroutine lgs_rb_y_2d(omega, ss, uu, ff, lo, ng)
    integer, intent(in) :: lo(:), ng
    real (kind = dp_t), intent(in) :: omega
    real (kind = dp_t), intent(in) :: ss(lo(1):,lo(2):,0:)
    real (kind = dp_t), intent(in) :: ff(lo(1):,lo(2):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t) :: wrk(size(ff,dim=2),4)
    integer i, j, hi(size(lo))

    hi = ubound(uu)-ng
    do i = lo(1), hi(1), 2
       do j = lo(2), hi(2)
          wrk(j,1:3) = ss(i,j,(/0,4,3/))
          wrk(j,4) = ff(i,j) -(ss(i,j,1)*uu(i+1,j) + ss(i,j,2)*uu(i-1,j))
       end do
       call dgtsl(wrk(:,2), wrk(:,1), wrk(:,3), wrk(:,4))
       do j = lo(2), hi(2)
          uu(i,j) = uu(i,j) + omega*(wrk(j,4)-uu(i,j))
       end do
    end do

    do i = lo(1)+1, hi(1), 2
       do j = lo(2), hi(2)
          wrk(j,1:3) = ss(i,j,(/0,4,3/))
          wrk(j,4) = ff(i,j) -(ss(i,j,1)*uu(i+1,j) + ss(i,j,2)*uu(i-1,j))
       end do
       call dgtsl(wrk(:,2), wrk(:,1), wrk(:,3), wrk(:,4))
       do j = lo(2), hi(2)
          uu(i,j) = uu(i,j) + omega*(wrk(j,4)-uu(i,j))
       end do
    end do

  end subroutine lgs_rb_y_2d

  subroutine dgtsl(c, d, e, b, lpivot)
    use bl_error_module
    real (kind = dp_t), intent(inout), dimension(:) :: c, d, e, b
    logical, optional, intent(in) :: LPIVOT

    !     dgtsl given a general tridiagonal matrix and a right hand
    !     side will find the solution.

    !     on entry

    !     n       integer
    !     is the order of the tridiagonal matrix.

    !     c       real (kind = dp_t)(n)
    !     is the subdiagonal of the tridiagonal matrix.
    !     c(2) through c(n) should contain the subdiagonal.
    !     on output c is destroyed.

    !     d       real (kind = dp_t)(n)
    !     is the diagonal of the tridiagonal matrix.
    !     on output d is destroyed.

    !     e       real (kind = dp_t)(n)
    !     is the superdiagonal of the tridiagonal matrix.
    !     e(1) through e(n-1) should contain the superdiagonal.
    !     on output e is destroyed.

    !     b       real (kind = dp_t)(n)
    !     is the right hand side vector.

    !     on return

    !     b       is the solution vector.

    !     linpack. this version dated 08/14/78 .
    !     jack dongarra, argonne national laboratory.

    !     no externals
    !     fortran dabs

    !     internal variables

    integer k, n
    real (kind = dp_t) t
    real (kind = dp_t), parameter :: ZERO = 0.0_dp_t
    logical :: lpv

    lpv = .TRUE.; if ( present(lpivot) ) lpv = lpivot

    !     begin block permitting ...exits to 100

    n = size(b)
    if ( n .lt. 1 ) return
    if ( n .eq. 1 ) then
       b(1) = b(1)/d(1)
       return
    end if

    c(1) = d(1)
    d(1) = e(1)
    e(1) = ZERO
    e(n) = ZERO

    do k = 1, n-1

       ! find the largest of the two rows

       if (LPV .and. abs(c(k+1)) .ge. abs(c(k))) then

          ! interchange row

          t = c(k+1)
          c(k+1) = c(k)
          c(k) = t
          t = d(k+1)
          d(k+1) = d(k)
          d(k) = t
          t = e(k+1)
          e(k+1) = e(k)
          e(k) = t
          t = b(k+1)
          b(k+1) = b(k)
          b(k) = t
       end if

       ! zero elements

       if (c(k) .eq. ZERO) then
          call bl_error("error 0")
       end if
       t = -c(k+1)/c(k)
       c(k+1) = d(k+1) + t*d(k)
       d(k+1) = e(k+1) + t*e(k)
       e(k+1) = ZERO
       b(k+1) = b(k+1) + t*b(k)
    end do
    if (c(n) .eq. ZERO) then
       call bl_error("error 1")
    end  if

    ! back solve

    b(n) = b(n)/c(n)
    b(n-1) = (b(n-1) - d(n-1)*b(n))/c(n-1)
    do k = n-2, 1, -1
       b(k) = (b(k) - d(k)*b(k+1) - e(k)*b(k+2))/c(k)
    end do

  end subroutine dgtsl

end module mg_smoother_module
