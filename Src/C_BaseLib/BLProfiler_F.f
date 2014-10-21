c-----------------------------------------------------------------------
      subroutine bl_proffortfuncstart(str)
      character*(*) str
      integer NSTR
      parameter (NSTR = 128)
      integer istr(NSTR)
      call blstr2int(istr, NSTR, str)
      call bl_proffortfuncstart_cpp(istr, NSTR)
      end
c-----------------------------------------------------------------------
      subroutine bl_proffortfuncstop(str)
      character*(*) str
      integer NSTR
      parameter (NSTR = 128)
      integer istr(NSTR)
      call blstr2int(istr, NSTR, str)
      call bl_proffortfuncstop_cpp(istr, NSTR)
      end
