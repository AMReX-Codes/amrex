      subroutine dsortem(ib,ie,a,iperm,b,c,d,e,f,g,h)
c-----------------------------------------------------------------------
c
c                      Quickersort Subroutine
c                      **********************
c
c This is a subroutine for sorting a real array in ascending order. This
c is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
c in collected algorithms of the ACM.
c
c The method used is that of continually splitting the array into parts
c such that all elements of one part are less than all elements of the
c other, with a third part in the middle consisting of one element.  An
c element with value t is chosen arbitrarily (here we choose the middle
c element). i and j give the lower and upper limits of the segment being
c split.  After the split a value q will have been found such that 
c a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
c performs operations on the two segments (i,q-1) and (q+1,j) as follows
c The smaller segment is split and the position of the larger segment is
c stored in the lt and ut arrays.  If the segment to be split contains
c two or fewer elements, it is sorted and another segment is obtained
c from the lt and ut arrays.  When no more segments remain, the array
c is completely sorted.
c
c
c INPUT PARAMETERS:
c
c   ib,ie        start and end index of the array to be sorteda
c   a            array, a portion of which has to be sorted.
c   iperm        0 no other array is permuted.
c                1 array b is permuted according to array a
c                2 arrays b,c are permuted.
c                3 arrays b,c,d are permuted.
c                4 arrays b,c,d,e are permuted.
c                5 arrays b,c,d,e,f are permuted.
c                6 arrays b,c,d,e,f,g are permuted.
c                7 arrays b,c,d,e,f,g,h are permuted.
c               >7 no other array is permuted.
c
c   b,c,d,e,f,g,h  arrays to be permuted according to array a.
c
c OUTPUT PARAMETERS:
c
c    a      = the array, a portion of which has been sorted.
c
c    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
c
c NO EXTERNAL ROUTINES REQUIRED:
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
c
c The dimensions for lt and ut have to be at least log (base 2) n
c
      integer   lt(64),ut(64),i,j,k,m,p,q
c
c Initialize:
c
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
c
c If this segment has more than two elements  we split it
c
 10   if (j-i-1) 100,90,15
c
c p is the position of an arbitrary element in the segment we choose the
c middle element. Under certain circumstances it may be advantageous
c to choose p at random.
c
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
c
c Start at the beginning of the segment, search for k such that a(k)>t
c
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
c
c Such an element has now been found now search for a q such that a(q)<t
c starting at the end of the segment.
c
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
c
c a(q) has now been found. we interchange a(q) and a(k)
c
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
c
c Update q and search for another pair to interchange:
c
      q = q-1
      go to 20
 50   q = k-1
 60   continue
c
c The upwards search has now met the downwards search:
c
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
c
c The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
c store the position of the largest segment in lt and ut
c
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
c
c Update m and split the new smaller segment
c
 80   m = m+1
      go to 10
c
c We arrive here if the segment has  two elements we test to see if
c the segment is properly ordered if not, we perform an interchange
c
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
c
c If lt and ut contain more segments to be sorted repeat process:
c
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      return
      end
