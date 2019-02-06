**********************************************************************  
**********************************************************************  
      subroutine detko(alf,phi,del)
*
* is called to determine coefficients of taylor series for interpolation
*	of deflection angle due to cells
*
* input:
*        alf(8)  - vector with four angles alphax(1), ..., alphay(4)    
*        del     - distance between two adjacent rays 
*
* output:
*        phi(8)  - coefficients
*
* called in:
*           main
*
* use of subroutines:
*           none
*
* use of functions:
*           none
*					October 31, 1991; J. Wambsganss
*
      implicit none
      integer debug
      double precision alf(8),phi(8),fac(4),del,oodel
      common/test/debug
*
      oodel  =  1.0d0/del
      fac(1) =  0.250d0
      fac(2) =  0.125d0*oodel
      fac(3) =  0.250d0*oodel*oodel
      fac(4) =  0.375d0*oodel*oodel*oodel
*
      phi(1) = fac(1)*
     &       (+alf(1)       +alf(3)       +alf(5)       +alf(7)        )
      phi(2) = fac(1)*
     &       (       +alf(2)       +alf(4)       +alf(6)       +alf(8) )
      phi(3) = fac(2)*
     &       (+alf(1)-alf(2)-alf(3)-alf(4)-alf(5)+alf(6)+alf(7)+alf(8) )
      phi(4) = fac(2)*
     &       (+alf(1)+alf(2)+alf(3)-alf(4)-alf(5)-alf(6)-alf(7)+alf(8) )
      phi(5) = fac(3)*
     &       (       -alf(2)       +alf(4)       -alf(6)       +alf(8) )
      phi(6) = fac(3)*
     &       (+alf(1)       -alf(3)       +alf(5)       -alf(7)        )
      phi(7) = fac(4)*
     &       (-alf(1)-alf(2)+alf(3)-alf(4)+alf(5)+alf(6)-alf(7)+alf(8) )
      phi(8) = fac(4)*
     &       (+alf(1)-alf(2)+alf(3)+alf(4)-alf(5)+alf(6)-alf(7)-alf(8) )
*
      return
      end
