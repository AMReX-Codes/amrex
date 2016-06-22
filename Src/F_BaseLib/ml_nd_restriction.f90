module ml_nodal_restriction_module

  use bl_constants_module
  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module

  implicit none

  private

  public :: periodic_add_copy, ml_nodal_injection, ml_nodal_restriction, ml_nodal_restriction_wrapper

contains

  subroutine ml_nodal_restriction(crse, fine, mm_fine, mm_crse, ir, inject, zero_only)
    use nodal_restriction_module
    use impose_neumann_bcs_module
    type(multifab),  intent(inout)        :: crse
    type(multifab),  intent(inout)        :: fine
    type(imultifab), intent(in   )        :: mm_fine
    type(imultifab), intent(in   )        :: mm_crse
    integer,         intent(in)           :: ir(:)
    logical,         intent(in), optional :: inject
    logical,         intent(in), optional :: zero_only

    integer             :: i, n, rmode, dm, ng
    integer             :: lo (get_dim(fine)), hi (get_dim(fine))
    integer             :: vlo(get_dim(fine)), vhi(get_dim(fine)) 
    integer             :: loc(get_dim(fine)), lof(get_dim(fine))
    integer             :: lom_fine(get_dim(fine)), lom_crse(get_dim(fine))
    logical             :: linject, lzero_only
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    integer,    pointer :: mp_fine(:,:,:,:), mp_crse(:,:,:,:)
    type(layout)        :: lacfine, laf
    type(multifab)      :: cfine
    type(imultifab)     :: mm_cfine
    type(mfiter)        :: mfi
    type(box)           :: bx, vbx

    type(bl_prof_timer), save :: bpt

    if ( ncomp(crse) .ne. ncomp(fine) ) then
       call bl_error('ml_nodal_restriction: crse & fine must have same # of components')
    end if

    call build(bpt, "ml_nodal_restriction")

    linject    = .false. ; if ( present(inject   ) ) linject    = inject
    lzero_only = .false. ; if ( present(zero_only) ) lzero_only = zero_only

    laf = get_layout(fine)

    call layout_build_coarse(lacfine, laf, ir)
    call multifab_build(cfine, lacfine, nc = ncomp(crse), ng = 0, nodal = nodal_flags(crse))
    call copy(cfine, crse)

    dm = get_dim(fine)

    if ( .not. linject ) then

       !$omp parallel private(mfi,n,i,bx,lo,hi,loc,lom_fine,cp,mp_fine)
       call mfiter_build(mfi,cfine,.true.)
       do n = 1, ncomp(cfine)
          do while(next_tile(mfi,i))
             bx = get_tilebox(mfi)

             lo       = lwb(bx)
             hi       = upb(bx)
             loc      = lwb(get_pbox(cfine,   i))
             lom_fine = lwb(get_pbox(mm_fine, i))

             cp      => dataptr(cfine,   i, n, 1)
             mp_fine => dataptr(mm_fine, i, 1, 1) ! mask has only 1 component
             select case (dm)
             case (1)
                call nodal_zero_1d(cp(:,1,1,1), loc, mp_fine(:,1,1,1), lom_fine, lo, hi, ir)
             case (2)
                call nodal_zero_2d(cp(:,:,1,1), loc, mp_fine(:,:,1,1), lom_fine, lo, hi, ir)
             case (3)
                call nodal_zero_3d(cp(:,:,:,1), loc, mp_fine(:,:,:,1), lom_fine, lo, hi, ir)
             end select
          end do
       end do
       !$omp end parallel

       call copy(crse, cfine)
       call setval(cfine, ZERO)
    end if

    if ( .not. lzero_only ) then

       if (.not. linject) then
          !$omp parallel do private(i,fp,mp_fine,lof,lom_fine,ng,n)
          do i = 1, nfabs(fine)
             fp      => dataptr(fine,    i)
             mp_fine => dataptr(mm_fine, i)
             lof      = lwb(get_pbox(fine,    i))
             lom_fine = lwb(get_pbox(mm_fine, i))
             ng       = lom_fine(1) - lof(1)
             do n = 1, ncomp(fine)
                select case (dm)
                case (1)
                   call impose_neumann_bcs_1d(fp(:,1,1,n),mp_fine(:,1,1,1),lom_fine,ng)
                case (2)
                   call impose_neumann_bcs_2d(fp(:,:,1,n),mp_fine(:,:,1,1),lom_fine,ng)
                case (3)
                   call impose_neumann_bcs_3d(fp(:,:,:,n),mp_fine(:,:,:,1),lom_fine,ng)
                end select
             end do
          end do
          !$omp end parallel do
       end if

       rmode = 0
       call imultifab_build(mm_cfine, lacfine, nc = ncomp(mm_crse), ng = 0, nodal = nodal_flags(mm_crse))
       call copy(mm_cfine, mm_crse)

       !$omp parallel private(mfi,n,i,bx,vbx,lo,hi,vlo,vhi,lof,loc,lom_crse,lom_fine) &
       !$omp& private(cp,fp,mp_crse,mp_fine)
       call mfiter_build(mfi,cfine,.true.)
       do n = 1, ncomp(cfine)
          do while(next_tile(mfi,i))
             bx  = get_tilebox(mfi)
             vbx = get_ibox(cfine, i)

             lo       = lwb(bx)
             hi       = upb(bx)
             vlo      = lwb(vbx)
             vhi      = upb(vbx)
             lof      = lwb(get_pbox(fine,    i))
             loc      = lwb(get_pbox(cfine,   i))
             lom_crse = lwb(get_pbox(mm_cfine,i))
             lom_fine = lwb(get_pbox(mm_fine, i))

             cp      => dataptr(cfine,   i, n, 1)
             fp      => dataptr(fine,    i, n, 1)
             mp_crse => dataptr(mm_cfine,i, 1, 1) ! mask has only 1 component
             mp_fine => dataptr(mm_fine, i, 1, 1) ! mask has only 1 component
             select case (dm)
             case (1)
                call nodal_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, &
                     mp_fine(:,1,1,1), lom_fine, &
                     mp_crse(:,1,1,1), lom_crse, lo, hi, vlo, vhi, ir, linject, rmode)
             case (2)
                call nodal_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, &
                     mp_fine(:,:,1,1), lom_fine, &
                     mp_crse(:,:,1,1), lom_crse, lo, hi, vlo, vhi, ir, linject, rmode)
             case (3)
                call nodal_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, &
                     mp_fine(:,:,:,1), lom_fine, &
                     mp_crse(:,:,:,1), lom_crse, lo, hi, vlo, vhi, ir, linject, rmode)
             end select
          end do
       end do
       !$omp end parallel

       call imultifab_destroy(mm_cfine)

       if ( linject ) then
          call multifab_copy(crse, cfine)
       else
          call multifab_copy(crse, cfine, filter = ml_restrict_copy_sum)
          call periodic_add_copy(crse,cfine,synced=.false.)
       end if
    end if

    call multifab_destroy(cfine)

    call destroy(bpt)

  end subroutine ml_nodal_restriction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_nodal_injection(crse, fine, ir)

    use nodal_restriction_module
    use impose_neumann_bcs_module
    type(multifab),  intent(inout)        :: crse
    type(multifab),  intent(inout)        :: fine
    integer,         intent(in)           :: ir(:)

    integer             :: i, n, dm
    integer             :: lo (get_dim(fine)), hi (get_dim(fine))
    integer             :: loc(get_dim(fine)), lof(get_dim(fine))
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(layout)        :: lacfine, laf
    type(multifab)      :: cfine
    type(mfiter)        :: mfi
    type(box)           :: bx

    type(bl_prof_timer), save :: bpt

    if ( ncomp(crse) .ne. ncomp(fine) ) then
       call bl_error('ml_nodal_injection: crse & fine must have same # of components')
    end if

    call build(bpt, "ml_nodal_injection")

    laf = get_layout(fine)

    call layout_build_coarse(lacfine, laf, ir)
    call multifab_build(cfine, lacfine, nc = ncomp(crse), ng = 0, nodal = nodal_flags(crse))
    call copy(cfine, crse)

    dm = get_dim(fine)

    !$omp parallel private(mfi,n,i,bx,lo,hi,lof,loc) &
    !$omp& private(cp,fp)
    call mfiter_build(mfi,cfine,.true.)
    do n = 1, ncomp(cfine)
       do while(next_tile(mfi,i))
          bx  = get_tilebox(mfi)

          lo       = lwb(bx)
          hi       = upb(bx)
          lof      = lwb(get_pbox(fine,    i))
          loc      = lwb(get_pbox(cfine,   i))

          cp      => dataptr(cfine,   i, n, 1)
          fp      => dataptr(fine,    i, n, 1)
          select case (dm)
          case (1)
             call nodal_injection_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, &
                  lo, hi, ir)
          case (2)
             call nodal_injection_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, &
                  lo, hi, ir)
          case (3)
             call nodal_injection_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, &
                  lo, hi, ir)
          end select
       end do
    end do
    !$omp end parallel

    call multifab_copy(crse, cfine)

    call multifab_destroy(cfine)

    call destroy(bpt)

  end subroutine ml_nodal_injection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine periodic_add_copy(dst,src,synced)

    use bl_prof_module

    type(multifab), intent(inout) :: dst
    type(multifab), intent(in   ) :: src
    logical,        intent(in   ) :: synced
    !
    ! if ( synced )
    !
    !   This version assumes that src IS synced up on each edge to start with  - for example,
    !   if a node in grid A on the lo-x side has value a, and the equivalent node in grid B 
    !   but on the hi-x side has value b, and the equivalent node in grid C on the hi-x side
    !   has value b also, then the final value of the nodes in A,B and C will be (a+b)
    !
    ! else
    !
    !   This version assumes that src is NOT synced up on each edge to start with  - for example,
    !             if a node in grid A on the lo-x side has value a and 
    !   the equivalent node in grid B on the hi-x side has value b and
    !   the equivalent node in grid C on the hi-x side has value c ,
    !   then the final value of the node in each grid A, B and C will be (a+b+c)
    !
    type(multifab)      :: temp_dst
    type(box)           :: domain,bxi,bxj,bx_to,bx_from
    type(box)           :: domain_edge_src, domain_edge_dst
    type(layout)        :: dstla, srcla
    real(dp_t), pointer :: ap(:,:,:,:), bp(:,:,:,:)
    integer             :: i,j,ii,jj,idir,jdir,kdir,proc,dm,nc
    logical             :: nodal(get_dim(dst))
    integer             :: shift_vector(3),dims(3),dlen(3)
    logical             :: pmask(get_dim(dst))

    integer,    parameter :: TAG  = 1111

    type(box_intersector), pointer :: bisrc(:), bidst(:)

    real(kind=dp_t), dimension(:,:,:,:), allocatable :: pt

    type(bl_prof_timer), save :: bpt

    if ( ncomp(dst) .ne. ncomp(src) ) then
       call bl_error('periodic_add_copy: src & dst must have same # of components')
    end if

    if ( .not. nodal_q(dst) ) call bl_error('periodic_add_copy(): dst NOT nodal')
    if ( .not. nodal_q(src) ) call bl_error('periodic_add_copy(): src NOT nodal')

    pmask = get_pmask(get_layout(dst))

    if ( all(pmask .eqv. .false.) ) return

    nc         = ncomp(dst)
    dm         = get_dim(dst)
    dims       = 1
    nodal      = .true.
    domain     = box_nodalize(get_pd(get_layout(dst)),nodal)
    dlen(1:dm) = extent(domain)

    dstla = get_layout(dst)
    srcla = get_layout(src)

    call build(bpt, "periodic_add_copy")

    if ( synced ) call multifab_build(temp_dst,dstla,nc,0,nodal)

    do kdir = -1,1

       if ( dm < 3  .and. kdir /= 0                         ) cycle
       if ( dm == 3 .and. (.not. pmask(dm)) .and. kdir /= 0 ) cycle

       if ( dm == 3 ) shift_vector(3) = kdir * (dlen(dm) - 1)

       do jdir = -1,1

          if ( .not. pmask(2) .and. jdir /= 0 ) cycle

          shift_vector(2) = jdir * (dlen(2) - 1)

          do idir = -1,1

             if ( .not. pmask(1) .and. idir /= 0                            ) cycle
             if ( dm == 2 .and. (idir == 0 .and. jdir == 0)                 ) cycle
             if ( dm == 3 .and. (idir == 0 .and. jdir == 0 .and. kdir == 0) ) cycle

             shift_vector(1) = idir * (dlen(1) - 1)

             domain_edge_src = intersection(domain,shift(domain, shift_vector))
             domain_edge_dst = shift(domain_edge_src,-shift_vector)

             if ( synced ) call setval(temp_dst,ZERO)
             !
             ! Add values from domain_edge_src side to domain_edge_dst side
             !
             bidst => layout_get_box_intersector(dstla, domain_edge_dst, nodal_la=nodal)

             do jj = 1, size(bidst)
                j     =  bidst(jj)%i
                bxj   =  bidst(jj)%bx
                bisrc => layout_get_box_intersector(srcla, domain_edge_src, nodal_la=nodal)
                do ii = 1, size(bisrc)
                   i = bisrc(ii)%i
                   if ( remote(dst%la,j) .and. remote(src%la,i) ) cycle
                   bxi     = shift(bisrc(ii)%bx,-shift_vector)
                   bx_from = intersection(bxi,bxj)
                   if ( empty(bx_from) ) cycle
                   bx_to   = bx_from
                   bx_from = shift(bx_from,shift_vector)
                   if ( local(dst%la,j) .and. local(src%la,i) ) then
                      if ( synced ) then
                         ap => dataptr(temp_dst,local_index(temp_dst,j),bx_to)
                         bp => dataptr(src,     local_index(src,     i), bx_from)
                         call cpy_d(ap,bp) ! ap = bp
                      else
                         ap => dataptr(dst,local_index(dst,j),bx_to)
                         bp => dataptr(src,local_index(src,i),bx_from)
                         call cpy_d(ap,bp,filter=ml_restrict_copy_sum) ! ap = ap + bp
                      end if
                   else if ( local(src%la,i) ) then
                      !
                      ! We own src.
                      !
                      bp => dataptr(src,local_index(src,i),bx_from)
                      call parallel_send(bp, get_proc(get_layout(dst),j), TAG)
                   else
                      !
                      ! We own dst.
                      !
                      dims(1:dm) = extent(bx_from)
                      proc = get_proc(get_layout(src),i)
                      allocate(pt(dims(1),dims(2),dims(3),nc))
                      call parallel_recv(pt, proc, TAG)
                      if ( synced ) then
                         ap => dataptr(temp_dst,local_index(temp_dst,j),bx_to)
                         call cpy_d(ap,pt) ! ap = pt
                      else
                         ap => dataptr(dst,local_index(dst,j),bx_to)
                         call cpy_d(ap,pt,filter=ml_restrict_copy_sum) ! ap = ap + pt
                      end if
                      deallocate(pt)
                   end if
                end do

                deallocate(bisrc)
             end do

             deallocate(bidst)

             if ( synced ) call plus_plus(dst,temp_dst)
          end do
       end do
    end do

    call clear_box_hash_bin(dstla%lap)
    call clear_box_hash_bin(srcla%lap)

    if ( synced ) call multifab_destroy(temp_dst)

    call destroy(bpt)

  end subroutine periodic_add_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! a wrapper for nodal restriction
 ! need to build mm_fine, mm_crse
 
 subroutine ml_nodal_restriction_wrapper(mla,mf,the_bc_tower,bc_comp)

   use nodal_stencil_bc_module

   type(ml_layout), intent(inout) :: mla
   type(multifab) , intent(inout) :: mf(:) ! n-level array of nodal multifabs
   type(bc_tower) , intent(in   ) :: the_bc_tower
   integer        , intent(in   ) :: bc_comp

   ! local
   type(imultifab) :: mm(mla%nlevel)
   integer, allocatable :: face_type(:,:,:)

   integer :: n,nlevs,i

   type(box) :: pd,bx,nbx
   integer :: id,lo_grid,hi_grid,lo_dom,hi_dom

   logical :: nodal_flag(mla%dim)

   type(layout) :: la

   nlevs = mla%nlevel

   nodal_flag(:) = .true.

   ! build the nodal mask for each level
   do n=1,nlevs

      call imultifab_build(mm(n), mla%la(n), 1, 0, nodal_flag)

      ! extract the problem domain
      pd = layout_get_pd(mla%la(n))

      ! allocate and fill face_type for the level
       allocate(face_type(nfabs(mf(n)),mla%dim,2))
       face_type = BC_INT
       do id = 1,mla%dim
          lo_dom = lwb(pd,id)
          hi_dom = upb(pd,id)
          do i = 1,nfabs(mf(n))
             lo_grid =  lwb(get_box(mf(n),i),id)
             if (lo_grid == lo_dom) then
                face_type(i,id,1) = &
                     the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,id,1,bc_comp)
             end if
             hi_grid = upb(get_box(mf(n),i),id)
             if (hi_grid == hi_dom) then
                face_type(i,id,2) = &
                     the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,id,2,bc_comp)
             end if
          end do
       end do


       ! build the mask for the level
       ! start by initializing to zero
       call setval(mm(n),BC_INT)  

       !$omp parallel do private(i,bx,nbx,la)
       do i = 1, nfabs(mf(n))
          bx  = get_box(mf(n),i)
          nbx = box_nodalize(bx,nodal_flag)
          la = mla%la(n)
          call stencil_set_bc_nodal(mla%dim, bx, nbx, i, mm(n), face_type, pd, la)
       end do
       !$omp end parallel do

       deallocate(face_type)

   end do

   ! do nodal restriction
   do n=nlevs,2,-1
      call ml_nodal_restriction(mf(n-1),mf(n),mm(n),mm(n-1),mla%mba%rr(n-1,:),inject=.true.)
   end do

   do n=1,nlevs
      call imultifab_destroy(mm(n))
   end do

 end subroutine ml_nodal_restriction_wrapper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_restrict_copy_sum(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     !
     ! out = out + in 
     !
     do n = 1, size(out,4)
        do k = 1, size(out,3)
           do j = 1, size(out,2)
              do i = 1, size(out,1)
                 out(i,j,k,n) = out(i,j,k,n) + in(i,j,k,n)
              end do
           end do
        end do
     end do
  end subroutine ml_restrict_copy_sum

end module ml_nodal_restriction_module
