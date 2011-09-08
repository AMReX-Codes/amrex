      subroutine ktsol(n,ns,nv,a,b,x,ktilt,maxeq)
c-----------------------------------------------------------------------
c
c Solution of a system of linear equations by gaussian elimination with
c partial pivoting.  Several right hand side matrices and several
c variables are allowed.
c
c
c         NOTE: All input matrices must be in double precision
c
c
c INPUT/OUTPUT VARIABLES:
c
c   n                Number of equations
c   ns               Number of right hand side matrices
c   nv               Number of variables.
c   a(n*n*nv)        left hand side matrices versus columnwise.
c   b(n*ns*nv)       input right hand side matrices.
c   x(n*ns*nv)       solution matrices.
c   ktilt            indicator of singularity
c                      =  0  everything is ok.
c                      = -1 n.le.1
c                      =  k  a null pivot appeared at the kth iteration.
c   tol              used in test for null pivot. depends on machine
c                      precision and can also be set for the tolerance
c                      of an ill-defined kriging system.
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      real*8 x(maxeq),a(maxeq*maxeq),b(maxeq)
c
c Make sure there are equations to solve:
c
      if(n.le.1) then
            ktilt = -1
            return
      endif
c
c Initialization:
c
      tol   = 0.1e-10
      ktilt = 0
      ntn   = n*n
      nm1   = n-1
c
c Triangulation is done variable by variable:
c
      do iv=1,nv
c
c Indices of location in vectors a and b:
c
            nva = ntn*(iv-1)
            nvb = n*ns*(iv-1)
c
c Gaussian elimination with partial pivoting:
c
            do k=1,nm1
                  kp1 = k+1
c
c Indice of the diagonal element in the kth row:
c
                  kdiag = nva+(k-1)*n+k
c
c Find the pivot - interchange diagonal element/pivot:
c
                  npiv = kdiag
                  ipiv = k
                  i1   = kdiag
                  do i=kp1,n
                        i1 = i1+1
                        if(abs(a(i1)).gt.abs(a(npiv))) then
                              npiv = i1
                              ipiv = i
                        endif
                  end do
                  t        = a(npiv)
                  a(npiv)  = a(kdiag)
                  a(kdiag) = t
c
c Test for singularity:
c
                  if(abs(a(kdiag)).lt.tol) then
                        ktilt=k
                        return
                  endif
c
c Compute multipliers:
c
                  i1 = kdiag
                  do i=kp1,n
                        i1    = i1+1
                        a(i1) = -a(i1)/a(kdiag)
                  end do
c
c Interchange and eliminate column per column:
c
                  j1 = kdiag
                  j2 = npiv
                  do j=kp1,n
                        j1    = j1+n
                        j2    = j2+n
                        t     = a(j2)
                        a(j2) = a(j1)
                        a(j1) = t
                        i1    = j1
                        i2    = kdiag
                        do i=kp1,n
                              i1    = i1+1
                              i2    = i2+1
                              a(i1) = a(i1)+a(i2)*a(j1)
                        end do
                  end do
c
c Interchange and modify the ns right hand matrices:
c
                  i1 = nvb+ipiv
                  i2 = nvb+k
                  do i=1,ns
                        t     = b(i1)
                        b(i1) = b(i2)
                        b(i2) = t
                        j1    = i2
                        j2    = kdiag
                        do j=kp1,n
                              j1    = j1+1
                              j2    = j2+1
                              b(j1) = b(j1)+b(i2)*a(j2)
                        end do
                        i1 = i1+n
                        i2 = i2+n
                  end do
            end do
c
c Test for singularity for the last pivot:
c
            kdiag = ntn*iv
            if(abs(a(kdiag)).lt.tol) then
                  ktilt = n
                  return
            endif
      end do
c
c End of triangulation. Now, solve back variable per variable:
c
      do iv=1,nv
c
c Indices of location in vectors a and b:
c
            nva  = ntn*iv
            nvb1 = n*ns*(iv-1)+1
            nvb2 = n*ns*iv
c
c Back substitution with the ns right hand matrices:
c
            do il=1,ns
                  do k=1,nm1
                        nmk = n-k
c
c Indice of the diagonal element of the (n-k+1)th row and of
c the (n-k+1)th element of the left hand side.
c
                        kdiag = nva-(n+1)*(k-1)
                        kb    = nvb2-(il-1)*n-k+1
                        b(kb) = b(kb)/a(kdiag)
                        t     = -b(kb)
                        i1    = kb
                        i2    = kdiag
                        do i=1,nmk
                              i1    = i1-1
                              i2    = i2-1
                              b(i1) = b(i1)+a(i2)*t
                        end do
                  end do
                  kdiag = kdiag-n-1
                  kb    = kb-1
                  b(kb) = b(kb)/a(kdiag)
            end do
c
c End of back substitution:
c
      end do
c
c Restitution of the solution:
c
      itot = n*ns*nv
      do i=1,itot
            x(i) = b(i)
      end do
c
c Finished:
c
      return
      end
