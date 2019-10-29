Program pointtest
   USE Mesh, ONLY : TPoint
   IMPLICIT NONE

   TYPE(TPoint):: p1
   p1 = TPoint(13.5d0, 0.2d0)
   if (p1%x /= 13.5d0 .OR. p1%y /= 0.2d0) call exit(1)
END PROGRAM pointtest
