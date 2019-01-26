
      subroutine ext_src_add(lo,hi,&
                           old_state,old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3,&
                           new_state,new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3,&
                           src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                           add,add_l1,add_l2,add_l3,add_h1,add_h2,add_h3, &
                           problo,dx,time,z,dt,bh_mass_acc,smbh_switch, &
                           spc_switch)

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : NVAR

      implicit none

      integer ,intent(in   ) :: lo(3),hi(3)
      integer ,intent(in   ) :: old_state_l1,old_state_l2,old_state_l3,old_state_h1,old_state_h2,old_state_h3
      integer ,intent(in   ) :: new_state_l1,new_state_l2,new_state_l3,new_state_h1,new_state_h2,new_state_h3
      integer ,intent(in   ) :: src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer ,intent(in   ) :: add_l1,add_l2,add_l3,add_h1,add_h2,add_h3
      real(rt),intent(in   ) :: old_state(old_state_l1:old_state_h1,old_state_l2:old_state_h2, &
                                                 old_state_l3:old_state_h3,NVAR)
      real(rt),intent(in   ) :: new_state(new_state_l1:new_state_h1,new_state_l2:new_state_h2, &
                                                 new_state_l3:new_state_h3,NVAR)
      real(rt),intent(  out) :: src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
      real(rt),intent(  out) :: add(add_l1:add_h1,add_l2:add_h2,add_l3:add_h3,NVAR)
      real(rt),intent(in   ) :: problo(3),dx(3),time,z,dt
      real(rt),intent(  out) :: bh_mass_acc
      integer ,intent(in   ) :: smbh_switch,spc_switch

      bh_mass_acc = 0.0d0

      src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.0d0
      add(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.0d0

      end subroutine ext_src_add
