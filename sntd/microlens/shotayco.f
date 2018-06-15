**********************************************************************  
**********************************************************************  
      subroutine shotayco(phi,del,lensmax,lensdata,lno,xxxx,yyyy,agam
     &,bgam,sigmac,sigmas,lensused,pix,ipix,pixminx,pixmaxx,pixminy
     &,pixmaxy,pix1,ipix1,factor2,iz,raydist2,xsidum,etadum,xxxxx,yyyyy
     &,alphax,alphay,ixdum,iydum,jxdum,jydum,nlens)
*
*  shoot with taylor coefficients in a square with factor2*factor2
*  rays and count rays 
*
* input:                                                                
*        phi(8)  - coefficients for Taylor expansion                    
*        del     - distance between center of this square and border
*        xxxx    - x coordinate of lower left corner of square          
*        yyyy    - y coordinate of lower left corner of square          
*        raydist2- distance of rays in this level
*        xxxxx   - (dummy) array with x-coordiantes of rays in square
*        yyyyy   - (dummy) array with y-coordiantes of rays in square
*
* output:                                                               
*        xsidum  - array with xsi-coordiantes of all rays in square     
*        etadum  - array with eta-coordiantes of all rays in square     
*
* called in:
*           main
*
* use of subroutines:                                                   
*           raycount                                                    
*
* use of functions:                                                     
*           none                                                        
*					October 31, 1991; J. Wambsganss
*
      implicit none                                                     
      real            times(20)
      integer*8 iz
      integer debug,lensmax,lno(lensmax),lensused,i1,indum,ipix,ipix1
     &,ix2,iy2,pix(ipix,ipix),pix1(ipix1,ipix1),factor2
     &,ixdum(factor2*factor2),iydum(factor2*factor2)
     &,jxdum(factor2*factor2),jydum(factor2*factor2),nlens
      external raycount
*
      double precision phi(8),del,sigmac,sigmas,xx,yy
     &,xx2,yy2,xxxx,yyyy,xxxxdel,yyyydel,z1,z2,z3,z4,lensdata(3,lensmax)
     &,dx,dy,modist2,xsidum(factor2,factor2),etadum(factor2,factor2)  
     &,alphax(factor2,factor2),alphay(factor2,factor2),pixminx,pixmaxx
     &,pixminy,pixmaxy,oopixfac,raydist2,pixminxt,pixminyt,pixmaxyt
     &,xxxxx(factor2,factor2),yyyyy(factor2,factor2),agam,bgam,pixmaxxt
*
      common/test/debug                                                 
      common/jkwtime/times
*
* center coordinates of this square:
*
      xxxxdel = xxxx + del
      yyyydel = yyyy + del
*
*
*
*
*
      do ix2= 1,factor2
         do iy2= 1,factor2
*
*    determine the new coordinates for rays in a square:
*    (for which the deflection angle should be calculated 
*     with Taylor series)
*
            xx             = ( float(ix2)-0.5  )*raydist2  - del    
            yy             = ( float(iy2)-0.5  )*raydist2  - del    
            xxxxx(ix2,iy2) = xxxxdel + xx
            yyyyy(ix2,iy2) = yyyydel + yy
*
            xx2 = xx*xx
            yy2 = yy*yy
            z1  = 0.5d0*(xx2 - yy2)
            z2  = xx*yy
            z3  = 0.5d0*xx*(xx2/3.0d0 - yy2    )
            z4  = 0.5d0*yy*(xx2     - yy2/3.0d0)
*
            alphax(ix2,iy2) =    phi(1)            +xx*phi(3) +yy*phi(4)
     &                       +z1*phi(5) +z2*phi(6) +z3*phi(7) +z4*phi(8)
*
            alphay(ix2,iy2) =               phi(2) -yy*phi(3) +xx*phi(4)
     &                       -z2*phi(5) +z1*phi(6) -z4*phi(7) +z3*phi(8)
         enddo
      enddo
*
      iz = iz+factor2*factor2
*
*
*
C$DIR SCALAR
      do ix2= 1,factor2
C$DIR SCALAR
         do iy2= 1,factor2
            do i1 = 1,lensused
*
* use this lens: position of ray: (xxx + del + xx, yyy + del + yy)
*
               indum           =   lno(i1)
               dx              = (xxxxx(ix2,iy2) - lensdata(1,indum))
               dy              = (yyyyy(ix2,iy2) - lensdata(2,indum))
               modist2         = lensdata(3,indum) / (dx*dx +  dy*dy)
               alphax(ix2,iy2) = alphax(ix2,iy2)  + dx*modist2
               alphay(ix2,iy2) = alphay(ix2,iy2)  + dy*modist2  
            enddo
         enddo
      enddo
*
	 pixminxt = pixminx 
	 pixminyt = pixminy
	 pixmaxxt = pixmaxx
	 pixmaxyt = pixmaxy
	 oopixfac = float(IPIX)/(pixmaxyt-pixminyt)
*
         do iy2= 1,factor2
            do ix2= 1,factor2
               xsidum(ix2,iy2)
     &  	=(agam-sigmac)*xxxxx(ix2,iy2)-alphax(ix2,iy2)
               etadum(ix2,iy2)
     &		=(bgam-sigmac)*yyyyy(ix2,iy2)-aLphay(ix2,iy2)
	    enddo
	 enddo
*
       call raycount(pix,ipix,pix1,ipix1,xsidum,etadum,factor2,pixminxt
     &,pixmaxxt,pixminyt,pixmaxyt,oopixfac,ixdum,iydum,jxdum,jydum)
*
      return
      end
