SUBROUTINE MATHERRQQ( name, length, info, retcode)
  USE DFLIB
  INTEGER(2) length, retcode
  CHARACTER(length) name
  RECORD /MTH$E_INFO/ info
  PRINT *, "Entered MATHERRQQ"
  PRINT *, "Failing function is: ", name
  PRINT *, "Error type is: ", info.errcode
  Call BL_ABORT()
  ! call Fort_abort()
  END
