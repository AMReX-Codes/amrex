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
	CALL bl_pp_get_string_cpp(ierr, ipp, istr, NSTR, ostr, NSTR)
	IF ( ierr ) THEN
	    CALL blint2str(ostr, iostr, NSTR)
      END IF
      END
