      subroutine ksol(nright,neq,nsb,a,r,s,ising)
c-----------------------------------------------------------------------
c
c                Solution of a System of Linear Equations
c                ****************************************
c
c
c
c INPUT VARIABLES:
c
c   nright,nsb       number of columns in right hand side matrix.
c                      for KB2D: nright=1, nsb=1
c   neq              number of equations
c   a()              upper triangular left hand side matrix (stored 
c                      columnwise)
c   r()              right hand side matrix (stored columnwise)
c                      for kb2d, one column per variable
c
c
c
c OUTPUT VARIABLES:
c
c   s()              solution array, same dimension as  r  above.
c   ising            singularity indicator
c                      0,  no singularity problem
c                     -1,  neq .le. 1
c                      k,  a null pivot appeared at the kth iteration
c
c
c
c PROGRAM NOTES:
c
c   1. Requires the upper triangular left hand side matrix.
c   2. Pivots are on the diagonal.
c   3. Does not search for max. element for pivot.
c   4. Several right hand side matrices possible.
c   5. USE for ok and sk only, NOT for UK.
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      real*8   a(*),r(*),s(*)
c
c If there is only one equation then set ising and return:
c
      if(neq.le.1) then
            ising = -1
            return
      endif
c
c Initialize:
c
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
c
c Start triangulation:
c
      do k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            do iv=1,nright
                  nm1=nm*(iv-1)
                  ii=kk+nn*(iv-1)
                  piv=1./a(ii)
                  lp=0
                  do i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
                        end do
                        do llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
                        end do
                  end do
            end do
      end do
c
c Error checking - singular matrix:
c
      ijm=ij-nn*(nright-1)
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
c
c Finished triangulation, start solving back:
c
      do iv=1,nright
            nm1=nm*(iv-1)
            ij=ijm+nn*(iv-1)
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
      end do
c
c Finished solving back, return:
c
      return
      end
