************************************************************************
************************************************************************
      subroutine rands(ar,cr1,cr2,var)
*
* generates a random number between 0 and 1 by selecting digits 5 to 13 
*  in a square r of a number that is square-rooted over and over;
*					(from Bohdan)
*
*  input:
*      ar       = a seed number
*      cr1, cr2 = auxillary seed numbers
*
*  output:
*      ar       = a modified seed number
*      var      = a random variable with a uniform probability distri
*                    in the range 0 to 1
*
* called in:
*           setstar
*
* use of subroutines:
*           none
*
* use of functions:
*           none
*					October 31, 1991; J. Wambsganss
*
      implicit none
      integer ia13,ib13,ia39,ib39
      double precision cr1,cr2,ar,br,f13,f26,al,bl,var,a13,b13,a13m,b13m
     &  ,a39,b39,r1,r2
*
      cr1   = cr1+1
      cr2   = cr2+1
      if(cr1.gt.1019)then
	 cr1 = 1.0d0
      endif
      if(cr2.gt.1056)then
	 cr2 = 1.0d0
      endif
      ar    = dsqrt(ar+cr1*cr2*1.765321d-5)
      br    = ar+dsqrt(1.1d0+0.1/ar)
      f13   = 2.0d0**13
      f26   = 2.0d0**26
*
      al    = dsqrt(ar)
      bl    = dsqrt(br)
      a13   = al*f13
      b13   = bl*f13
      ia13  = idint(a13)
      ib13  = idint(b13)
      a13m  = a13-ia13
      b13m  = b13-ib13
      a39   = a13m*f26
      b39   = b13m*f26
      ia39  = idint(a39)
      ib39  = idint(b39)
      r1    = ia39/f26
      r2    = ib39/f26
*
      var   = r1+r2/f26
      ar    = r2+r1/f26
*
      return
      end
