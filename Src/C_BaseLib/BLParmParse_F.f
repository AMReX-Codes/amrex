c
c $Id: BLParmParse_F.f,v 1.7 2003-01-03 19:25:11 car Exp $
c
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_new(ipp, str)
      INTEGER ipp
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_new_cpp(ipp, istr, NSTR)
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_record_new(ipp, ippr, str)
      INTEGER ipp, ippr
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_record_new_cpp(ipp, ippr, istr, NSTR)
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_get_int(ierr, ipp, str, ival)
      INTEGER ierr, ipp
      INTEGER ival
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
C
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_get_int_cpp(ierr, ipp, istr, NSTR, ival)
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_get_int_n(ierr, ipp, str, ival, nval)
      INTEGER ierr, ipp, nval
      INTEGER ival(*)
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
C
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_get_int_n_cpp(ierr, ipp, istr, NSTR, ival, nval)
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_get_logical(ierr, ipp, str, lval)
      INTEGER ierr, ipp
      LOGICAL lval
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
      INTEGER ival
C
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_get_logical_cpp(ierr, ipp, istr, NSTR, ival)
      IF ( IERR .NE. 0 ) THEN
          IF ( ival .NE. 0 ) THEN
              lval = .TRUE.
          ELSE
              lval = .FALSE.
          END IF
      END IF
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_get_real(ierr, ipp, str, rval)
      INTEGER ierr, ipp
      REAL rval
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
C
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_get_real_cpp(ierr, ipp, istr, NSTR, rval)
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_get_real_n(ierr, ipp, str, rval,nval)
      INTEGER ierr, ipp, nval
      REAL rval(*)
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
C
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_get_real_n_cpp(ierr, ipp, istr, NSTR, rval,nval)
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_get_double(ierr, ipp, str, dval)
      INTEGER ierr, ipp
      DOUBLE PRECISION dval
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
C
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_get_double_cpp(ierr, ipp, istr, NSTR, dval)
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_get_double_n(ierr, ipp, str, dval, nval)
      INTEGER ierr, ipp, nval
      DOUBLE PRECISION dval
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
C
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_get_double_n_cpp(ierr, ipp, istr, NSTR, dval,nval)
      END
c-----------------------------------------------------------------------
      SUBROUTINE bl_pp_get_string(ierr, ipp, str, ostr)
      INTEGER ierr, ipp
      CHARACTER*(*) ostr
      CHARACTER*(*) str
      INTEGER NSTR
      PARAMETER (NSTR = 128)
      INTEGER istr(NSTR)
      INTEGER iostr(NSTR)
C
      CALL blstr2int(istr, NSTR, str)
      CALL bl_pp_get_string_cpp(ierr, ipp, istr, NSTR, iostr, NSTR)
      IF ( ierr .ne. 0 ) THEN
          CALL blint2str(ostr, iostr, NSTR)
      END IF
      END
