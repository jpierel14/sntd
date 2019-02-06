**********************************************************************  
**********************************************************************  
      subroutine raycount(pix,ipix,pix1,ipix1,xsidum,etadum,factor2
     &,pixminx,pixmaxx,pixminy,pixmaxy,oopixfac,ix,iy,jx,jy)
*
* counts rays in pixels 
*
* input :
*        pix(ipix,ipix)               ---> pixel matrix
*        pix1(ipix,ipix)              ---> pixel matrix
*        xsi,eta                      ---> ray coordinates
*
* output: pix(ipix,ipix)              ---> pixel matrix
*         pix1(ipix,ipix)             ---> pixel matrix
*
* called in:
*           shotayco
*
* use of subroutines:
*           none
*
* use of functions:
*           none
*					October 31, 1991; J. Wambsganss
*
      implicit none
      integer factor2,debug,ix(factor2*factor2),ipix,ipix1,iz,iadd
     &,iy(factor2*factor2),jx(factor2*factor2),jy(factor2*factor2)
     &,kkk,ixdum,iydum
      integer*4 pix(ipix,ipix),pix1(ipix1,ipix1)
      double precision oopixfac,pixminx,pixmaxx,pixminy,pixmaxy,pixdifh
     &,xsidum(factor2*factor2),etadum(factor2*factor2),facipix
     &,dummy1,dummy2,dummy3,dummy4,xsidumdum,etadumdum
*
      common/test/debug
      common/pix/facipix,pixdifh
*
      dummy1 =  -pixminx*oopixfac + 0.5d0
      dummy2 =  -pixminy*oopixfac + 0.5d0
      dummy3 = (-pixminx+pixdifh)*facipix + 0.5d0
      dummy4 = (-pixminy+pixdifh)*facipix + 0.5d0
*
* "dirty programming": one loop instead of two because of faster 
*         vectorization
*
      kkk = 0
      do iz= 1,factor2*factor2
*
* determination of number of pixel hit by this ray (in array pix):  
*   (plus rays outside receiving field in rows/columns 0,ipix+1)
*
         xsidumdum  = xsidum(iz)*oopixfac+dummy1
	 ixdum      = nint(xsidumdum)
         etadumdum  = etadum(iz)*oopixfac+dummy2
         iydum      = nint(etadumdum)
*
         if(ixdum.lt.1.or.ixdum.gt.ipix.or
     &	   .iydum.lt.1.or.iydum.gt.ipix)then
	    goto 150
         else
	    kkk = kkk + 1
            ix(kkk)     = ixdum
            iy(kkk)     = iydum 
	 endif
  150    continue
*
      enddo
*
* now counting (with non-vectorizing loop)
*
      do iz=1,kkk
         iadd = 1
         pix(ix(iz),iy(iz)) = pix(ix(iz),iy(iz)) + iadd
      enddo
*
*
      do iz= 1,factor2*factor2
*
* for receiving field twice the sidelength:  (control!)
* determination of number of pixel hit by this ray (in array pix1): 
* (plus rays outside receiving field in rows/columns 0,ipix1+1)
*
         jx(iz)  = nint(xsidum(iz)*facipix + dummy3)
         jy(iz)  = nint(etadum(iz)*facipix + dummy4)
      enddo
*
      do iz=1,factor2*factor2
	 ixdum = jx(iz)
	 iydum = jy(iz)
	 if(ixdum.lt.1.or.ixdum.gt.ipix1.or
     &	   .iydum.lt.1.or.iydum.gt.ipix1)then
	    goto 250
	 else
            pix1(ixdum,iydum) = pix1(ixdum,iydum) + 1
	 endif
  250    continue
      enddo
*
      return
      end
