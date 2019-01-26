*DECK XSETUN
      SUBROUTINE XSETUN (LUN)
C***BEGIN PROLOGUE  XSETUN
C***PURPOSE  Reset the logical unit number for error messages.
C***CATEGORY  R3B
C***TYPE      ALL (XSETUN-A)
C***KEYWORDS  ERROR CONTROL
C***DESCRIPTION
C
C   XSETUN sets the logical unit number for error messages to LUN.
C
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETUN
C
C Subroutines called by XSETUN.. None
C Function routine called by XSETUN.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER LUN, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETUN
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETUN ----------------------
      END
