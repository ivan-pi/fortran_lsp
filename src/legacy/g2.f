      SUBROUTINE G2    (CTERM,STERM,X,Y)
c
C  APPLY THE ROTATION COMPUTED BY G1 TO (X,Y).
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1972 DEC 15, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
c     ------------------------------------------------------------------
      double precision CTERM, STERM, X, XR, Y
c     ------------------------------------------------------------------
      XR=CTERM*X+STERM*Y
      Y=-STERM*X+CTERM*Y
      X=XR
      RETURN
      END