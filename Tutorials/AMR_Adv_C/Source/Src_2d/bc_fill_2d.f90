   ! ::: -----------------------------------------------------------
 
   subroutine xvelfill(xvel,xvel_l1,xvel_l2,xvel_h1,xvel_h2,&
                       domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: xvel_l1,xvel_l2,xvel_h1,xvel_h2
     integer :: bc(2,2,*)
     integer :: domlo(2), domhi(2)
     double precision dx(2), xlo(2), time
     double precision xvel(xvel_l1:xvel_h1,xvel_l2:xvel_h2)
 
     call filcc(xvel,xvel_l1,xvel_l2,xvel_h1,xvel_h2,domlo,domhi,dx,xlo,bc)
 
   end subroutine xvelfill

   ! ::: -----------------------------------------------------------
 
   subroutine yvelfill(yvel,yvel_l1,yvel_l2,yvel_h1,yvel_h2,&
                       domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: yvel_l1,yvel_l2,yvel_h1,yvel_h2
     integer :: bc(2,2,*)
     integer :: domlo(2), domhi(2)
     double precision dx(2), xlo(2), time
     double precision yvel(yvel_l1:yvel_h1,yvel_l2:yvel_h2)
 
     call filcc(yvel,yvel_l1,yvel_l2,yvel_h1,yvel_h2,domlo,domhi,dx,xlo,bc)
 
   end subroutine yvelfill

   ! ::: -----------------------------------------------------------
 
   subroutine denfill(den,den_l1,den_l2,den_h1,den_h2,&
                      domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: den_l1,den_l2,den_h1,den_h2
     integer :: bc(2,2,*)
     integer :: domlo(2), domhi(2)
     double precision dx(2), xlo(2), time
     double precision den(den_l1:den_h1,den_l2:den_h2)
 
     call filcc(den,den_l1,den_l2,den_h1,den_h2,domlo,domhi,dx,xlo,bc)
 
   end subroutine denfill

   ! ::: -----------------------------------------------------------
 
   subroutine advfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,&
                       domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: adv_l1,adv_l2,adv_h1,adv_h2
     integer :: bc(2,2,*)
     integer :: domlo(2), domhi(2)
     double precision dx(2), xlo(2), time
     double precision adv(adv_l1:adv_h1,adv_l2:adv_h2)
 
     call filcc(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,dx,xlo,bc)
 
   end subroutine advfill

