c-----------------------------------------------------------------------
      subroutine bl_error(str)
      character*(*) str
      integer NSTR
      parameter (NSTR = 128)
      integer istr(NSTR)
      call blstr2int(istr, NSTR, str)
      call bl_error_cpp(istr, NSTR)
      end
c-----------------------------------------------------------------------
      subroutine bl_warning(str)
      character*(*) str
      integer NSTR
      parameter (NSTR = 128)
      integer istr(NSTR)
      call blstr2int(istr, NSTR, str)
      call bl_warning_cpp(istr, NSTR)
      end
c-----------------------------------------------------------------------
      subroutine bl_abort(str)
      character*(*) str
      integer NSTR
      parameter (NSTR = 128)
      integer istr(NSTR)
      call blstr2int(istr, NSTR, str)
      call bl_abort_cpp(istr, NSTR)
      end
