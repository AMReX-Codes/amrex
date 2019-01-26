      ! Compute the external source terms. This is called twice in the evolution:
      !
      ! 1) For the predictor, it is called with (old, old) states
      !    This is also used in the first pass of the conservative update
      !    (adding dt * S there)
      !
      ! 2) Next we correct the source terms in the conservative update to
      ! time-center them.  Here we call ext_src(old, new), and then
      ! in time_center_source_terms we subtract off 1/2 of the first S
      ! and add 1/2 of the new S.
      !
      ! Therefore, to get a properly time-centered source, generally
      ! speaking, you always want to use the "new" state here.  That
      ! will be the time n state in the first call and the n+1 in the
      ! second call.
      !

      subroutine ext_src(lo,hi,&
                         old_state,old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3,&
                         new_state,new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3,&
                         old_diag ,old_diag_l1,old_diag_l2,old_diag_l3,old_diag_h1,old_diag_h2,old_diag_h3,&
                         new_diag ,new_diag_l1,new_diag_l2,new_diag_l3,new_diag_h1,new_diag_h2,new_diag_h3,&
                         src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,problo,dx,time,z,dt) &
                         bind(C, name="fort_ext_src")

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : NVAR
      use turbforce_module  , only : force_scale

      implicit none

      integer ,intent(in   ) :: lo(3),hi(3)
      integer ,intent(in   ) :: old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3
      integer ,intent(in   ) :: new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3
      integer ,intent(in   ) :: old_diag_l1,old_diag_l2,old_diag_l3,old_diag_h1,old_diag_h2,old_diag_h3
      integer ,intent(in   ) :: new_diag_l1,new_diag_l2,new_diag_l3,new_diag_h1,new_diag_h2,new_diag_h3
      integer ,intent(in   ) :: src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      real(rt),intent(in   ) :: old_state(old_state_l1:old_state_h1,old_state_l2:old_state_h2, &
                                                 old_state_l3:old_state_h3,NVAR)
      real(rt),intent(in   ) :: new_state(new_state_l1:new_state_h1,new_state_l2:new_state_h2, &
                                                 new_state_l3:new_state_h3,NVAR)
      real(rt),intent(in   ) :: old_diag(old_diag_l1:old_diag_h1,old_diag_l2:old_diag_h2, &
                                                 old_diag_l3:old_diag_h3,NVAR)
      real(rt),intent(in   ) :: new_diag(new_diag_l1:new_diag_h1,new_diag_l2:new_diag_h2, &
                                                 new_diag_l3:new_diag_h3,NVAR)
      real(rt),intent(  out) :: src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
      real(rt),intent(in   ) :: problo(3),dx(3),time,z,dt

      if (force_scale > 0.0d0) then
          call ext_src_force(lo,hi,&
                             old_state,old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3,&
                             new_state,new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3,&
                             old_diag ,old_diag_l1,old_diag_l2,old_diag_l3,old_diag_h1,old_diag_h2,old_diag_h3,&
                             new_diag ,new_diag_l1,new_diag_l2,new_diag_l3,new_diag_h1,new_diag_h2,new_diag_h3,&
                             src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,problo,dx,time,z,dt)
      else
          src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0
      end if
 
      end subroutine ext_src

