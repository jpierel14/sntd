**********************************************************************  
**********************************************************************  
      subroutine numstar(nlens,lensmax,lensdata,levelmax,indlens)
*
* for every lens: calculates number of cell in every level and 
* compresses these numbers in one integer in vector indlens(nlens),
* then sorts lenses according to this number
*
* (upper left - 0, upper right - 1, lower left - 2, lower right - 3)
*
*
*    example:
*
*         level:       1      2      3       4      5     ..    levelmax
*    position  :       1      0      2       2      1     ..       2    
*  indlens(.) =        1*4**(levelmax-1) + 
*			      0*4**(levelmax-2)+
*			     	      2*4**(levelmax-3)+
*   							... 
*   							        +2*4**0
* input :
*        nlens           - number of lenses
*        lensmax         - maximum possible number of lenses
*        lensdata(1:3,i) - lens coordinates
*        levelmax   	 - maximum possible number of levels
*
*
* output: indlens(1:nlens) - numbers of cells in different
*                                levels, compressed (cf. example above)
*
*
* called in:
*           main
*
* use of subroutines:
*           sort
*
* use of functions:
*           none
*					October 31, 1991; J. Wambsganss
*
      implicit none
      external sort
      integer debug,lensmax,nlens,ilens,indlens(lensmax),index,ndum
     &,levelmax,fourpow(30),i,ii
*
      double precision xl,yl,lensdata(3,lensmax),stamin,stadif
     &,dumxmin,dumymin,dumdif,xratio, yratio,xcomp,ycomp
*
      common/leng/stamin,stadif
      common/test/debug
*
*
*
      xcomp   = 0.5
      ycomp   = 0.5
*
      do i = 1,levelmax  
        fourpow(i) = 4**(i-1)
      enddo
*
* loop over all lenses: numbers are calculated according to the cell
* position in different levels:
*
      do ilens = 1,nlens
	 xl      = lensdata(1,ilens)
	 yl      = lensdata(2,ilens)
  	 index   = 0
      	 dumxmin = stamin
      	 dumymin = stamin
      	 dumdif  = stadif
*
         do ii = 1,levelmax
            xratio = (xl-dumxmin)/dumdif
            yratio = (yl-dumymin)/dumdif
            dumdif = dumdif/2.0
            ndum   = 0
            if(xratio.gt.xcomp)then
                ndum    = ndum + 1
                dumxmin = dumxmin + dumdif
            endif
            if(yratio.lt.ycomp)then
                ndum    = ndum + 2
            else
                dumymin = dumymin + dumdif
            endif
            index = index + ndum*fourpow(levelmax-ii+1)
	 enddo
*
 	 indlens(ilens) = index
      enddo
*
*  sorts all lenses accoring to indlens(i): 
*
      if(nlens.gt.1)then
 	 call sort(nlens,indlens,lensdata,3)
      endif
*
      return
      end
