#include "FLUSH_F.H"

c-----------------------------------------------------------------------
      subroutine hgdebuginit(procno, numpc)
      integer procno, numpc
#include "hgdebug_F.H"
      character*24 hgdump
      dbunit = procno+10
      dbnorm = .FALSE.
      dbmat  = .FALSE.
#if   !defined(HG_DEBUG)
      dbon   = .FALSE.
#else
      dbon   = .TRUE.
      write(unit=hgdump,fmt='(a,i3.3,''_'',i3.3)')
     &     'fgu', numpc, procno
      open(unit=dbunit, file=hgdump, status='replace')
#endif
      end
c-----------------------------------------------------------------------
      subroutine hgdebugon
#include "hgdebug_F.H"
      dbon   = .TRUE.
      end
c-----------------------------------------------------------------------
      subroutine hgdebugoff
#include "hgdebug_F.H"
      dbon   = .FALSE.
      end
c-----------------------------------------------------------------------
      subroutine hgdebugctl(on, norm, mat)
      integer norm, mat, on
#include "hgdebug_F.H"
      dbon   = .TRUE.
      if ( on   .eq. 0 ) dbon   = .FALSE.
      dbnorm = .TRUE.
      if ( norm .eq. 0 ) dbnorm = .FALSE.
      dbmat  = .TRUE.
      if ( mat  .eq. 0 ) dbmat  = .FALSE.
      end
c-----------------------------------------------------------------------
      subroutine hgdebuglog(msg)
      character*(*) msg
#include "hgdebug_F.H"
      if ( .NOT. dbon ) return
      write(unit=dbunit, fmt=*) msg
      end
c-----------------------------------------------------------------------
      subroutine hgdebugs(name, s)
      character*(*) name, s
#include "hgdebug_F.H"
      if ( .NOT. dbon ) return
      write(unit=dbunit, fmt=*) name, s
      end
c-----------------------------------------------------------------------
      subroutine hgdebugi(name, i)
      integer i
      character*(*) name
#include "hgdebug_F.H"
      if ( .NOT. dbon ) return
      write(unit=dbunit, fmt=*) name, i
      call FORT_FLUSH(dbunit)
      end
c-----------------------------------------------------------------------
      subroutine hgdebugd(name, d)
      double precision d
      character*(*) name
#include "hgdebug_F.H"
      if ( .NOT. dbon ) return
      write(unit=dbunit, fmt=*) name, d
      call FORT_FLUSH(dbunit)
      end
#if (BL_SPACEDIM == 2)
c-----------------------------------------------------------------------
      subroutine hgdebugmat(name,
     &     a,   al0,al1,ah0,ah1,
     &     il0,il1,ih0,ih1, nc)
      character*(*) name
      integer al0,al1,ah0,ah1,nc
      integer il0,il1,ih0,ih1
      double precision a(al0:ah0,al1:ah1,1:nc)
      double precision norm(4)
      integer i,j,n
#include "hgdebug_F.H"
      if ( .NOT. dbon ) return
      if ( dbnorm ) then
         do i = 1, 4
            norm(i) = 0.0D0
         end do
         do n = 1, nc
            do i = il0, ih0
               do j = il1, ih1
                  norm(1) = norm(1) + a(i,j,n)
                  norm(3) = norm(3) + a(i,j,n)**2
               end do
            end do
         end do
         norm(3) = sqrt(norm(3))
         do n = 1, nc
            do i = al0, ah0
               do j = al1, ah1
                  norm(2) = norm(2) + a(i,j,n)
                  norm(4) = norm(4) + a(i,j,n)**2
               end do
            end do
         end do
         norm(4) = sqrt(norm(4))
      end if
      write(unit=dbunit, fmt=*) name, '(',
     &     il0, ':', ih0, ',',
     &     il1, ':', ih1, ',',
     &     '1', ':', nc,
     &     ')'
      if ( dbnorm ) then
         write(unit=dbunit, fmt=*) name,
     &        'n1(', norm(1), '/', norm(3), ')',
     &        'n2(', norm(2), '/', norm(4), ')'
      end if
      if ( dbmat ) then
         do n = 1, nc
            do i = il0, ih0
               do j = il1, ih1
                  write(unit=dbunit, fmt=*) name,
     &                 '(', i,j,n,')', a(i,j,n)
               end do
            end do
         end do
      end if
      call FORT_FLUSH(dbunit)
      end
#else
c-----------------------------------------------------------------------
      subroutine hgdebugmat(name,
     &     a,   al0,al1,al2,ah0,ah1,ah2,
     &          il0,il1,il2,ih0,ih1,ih2, nc)
      character*(*) name
      integer al0,al1,al2,ah0,ah1,ah2,nc
      integer il0,il1,il2,ih0,ih1,ih2
      double precision a(al0:ah0,al1:ah1,al2:ah2,1:nc)
      double precision norm(4)
      integer i,j,k,n
#include "hgdebug_F.H"
      if ( .NOT. dbon ) return
      if ( dbnorm ) then
         do i = 1, 4
            norm(i) = 0.0D0
         end do
         do n = 1, nc
            do i = il0, ih0
               do j = il1, ih1
                  do k = il2, ih2
                     norm(1) = norm(1) + a(i,j,k,n)
                     norm(3) = norm(3) + a(i,j,k,n)**2
                  end do
               end do
            end do
         end do
         norm(3) = sqrt(norm(3))
         do n = 1, nc
            do i = al0, ah0
               do j = al1, ah1
                  do k = al2, ah2
                     norm(2) = norm(2) + a(i,j,k,n)
                     norm(4) = norm(4) + a(i,j,k,n)**2
                  end do
               end do
            end do
         end do
         norm(4) = sqrt(norm(4))
      end if
      write(unit=dbunit, fmt=*) name, ': P(',
     &     al0, ':', ah0, ',',
     &     al1, ':', ah1, ',',
     &     al2, ':', ah2, ',',
     &     '1', ':', nc,
     &     ')'
      write(unit=dbunit, fmt=*) name, ': R(',
     &     il0, ':', ih0, ',',
     &     il1, ':', ih1, ',',
     &     il2, ':', ih2, ',',
     &     ')'
      if ( dbnorm ) then
         write(unit=dbunit, fmt=*) name,
     &        'n1(', norm(1), '/', norm(3), ')',
     &        'n2(', norm(2), '/', norm(4), ')'
      end if
      if ( dbmat ) then
         do n = 1, nc
            do i = il0, ih0
               do j = il1, ih1
                  do k = il2, ih2
                     write(unit=dbunit, fmt=*) name,
     &                    '(', i,j,k,n,')', a(i,j,k,n)
                  end do
               end do
            end do
         end do
      end if
      call FORT_FLUSH(dbunit)
      end
#endif
