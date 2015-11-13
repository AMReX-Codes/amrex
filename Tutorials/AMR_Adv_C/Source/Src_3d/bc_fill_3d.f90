   ! ::: -----------------------------------------------------------
 
   subroutine xvelfill(xvel,xvel_l1,xvel_l2,xvel_l3,xvel_h1,xvel_h2,xvel_h3,&
                       domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: xvel_l1,xvel_l2,xvel_l3,xvel_h1,xvel_h2,xvel_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision dx(3), xlo(3), time
     double precision xvel(xvel_l1:xvel_h1,xvel_l2:xvel_h2,xvel_l3:xvel_h3)
 
     call filcc(xvel,xvel_l1,xvel_l2,xvel_l3,xvel_h1,xvel_h2,xvel_h3,domlo,domhi,dx,xlo,bc)
 
   end subroutine xvelfill

   ! ::: -----------------------------------------------------------
 
   subroutine yvelfill(yvel,yvel_l1,yvel_l2,yvel_l3,yvel_h1,yvel_h2,yvel_h3,&
                       domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: yvel_l1,yvel_l2,yvel_l3,yvel_h1,yvel_h2,yvel_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision dx(3), xlo(3), time
     double precision yvel(yvel_l1:yvel_h1,yvel_l2:yvel_h2,yvel_l3:yvel_h3)
 
     call filcc(yvel,yvel_l1,yvel_l2,yvel_l3,yvel_h1,yvel_h2,yvel_h3,domlo,domhi,dx,xlo,bc)
 
   end subroutine yvelfill

   ! ::: -----------------------------------------------------------
 
   subroutine zvelfill(zvel,zvel_l1,zvel_l2,zvel_l3,zvel_h1,zvel_h2,zvel_h3,&
                       domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: zvel_l1,zvel_l2,zvel_l3,zvel_h1,zvel_h2,zvel_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision dx(3), xlo(3), time
     double precision zvel(zvel_l1:zvel_h1,zvel_l2:zvel_h2,zvel_l3:zvel_h3)
 
     call filcc(zvel,zvel_l1,zvel_l2,zvel_l3,zvel_h1,zvel_h2,zvel_h3,domlo,domhi,dx,xlo,bc)

   end subroutine zvelfill

   ! ::: -----------------------------------------------------------
 
   subroutine denfill(den,den_l1,den_l2,den_l3,den_h1,den_h2,den_h3,&
                       domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: den_l1,den_l2,den_l3,den_h1,den_h2,den_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision dx(3), xlo(3), time
     double precision den(den_l1:den_h1,den_l2:den_h2,den_l3:den_h3)
 
     call filcc(den,den_l1,den_l2,den_l3,den_h1,den_h2,den_h3,domlo,domhi,dx,xlo,bc)

   end subroutine denfill

   ! ::: -----------------------------------------------------------
 
   subroutine advfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,&
                       domlo,domhi,dx,xlo,time,bc)
 
     use probdata_module
     implicit none
     include 'bc_types.fi'
 
     integer :: adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
     integer :: bc(3,2,*)
     integer :: domlo(3), domhi(3)
     double precision dx(3), xlo(3), time
     double precision adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
 
     call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3,domlo,domhi,dx,xlo,bc)

   end subroutine advfill
