*DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH ()
C***BEGIN PROLOGUE  DUMACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        DOUBLE PRECISION  A, DUMACH
C        A = DUMACH()
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   930216  DATE WRITTEN
C   930818  Added SLATEC-format prologue.  (FNF)
C***END PROLOGUE  DUMACH
C
C*Internal Notes:
C-----------------------------------------------------------------------
C Subroutines/functions called by DUMACH.. None
C-----------------------------------------------------------------------
C**End
C
      DOUBLE PRECISION U
      DUMACH = EPSILON(U)
C----------------------- End of Function DUMACH ------------------------
      END
