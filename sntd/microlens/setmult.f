**********************************************************************  
**********************************************************************  
      subroutine setmult(lensdata,lensmax,nlens,cell,cellmax,ncell,i1   
     &	,i2)
*
* calculates quadrupole and higher moments for every cell that contains 
* more than one lens;  called once, after the set-up of all lenses
*
* input :
*         lensdata(3,.)            data of lenses
*         lensmax                  maximal number of lenses
*         ncell                  actual number of this cell
*         cell(14,.)                data of cells                       
*         cellmax                  maximal number of cells
*         nlens                    total number of lenses
*         i1                       position of first lens in this cell
*         i2                       position of last lens in this cell
*
* output: cell(5,k),cell(6,k)   quadrupole terms of cell #k             
*         cell(7,k),cell(8,k)   octopole   terms of cell #k             
*         cell(9,k),cell(10,k)   16-pole    terms of cell #k            
*         cell(11,k),cell(12,k)   32-pole  terms of cell #k             
*         cell(13,k),cell(14,k)   64-pole  terms of cell #k             
*
* called in:                                                            
*           setcell                                                     
*
* use of subroutines:                                                   
*           none
*
* use of functions:                                                     
*           none
*					October 31, 1991; J. Wambsganss
*
      implicit none                                                     
      integer cellmax,lensmax,nlens,ncell,debug,i1,i2,ii
      double precision cell(14,cellmax),lensdata(3,lensmax),xdum,xdum2
     &,xdum3,xdum4,xdum5,xdum6,ydum,ydum2,ydum3,ydum4,ydum5,ydum6
      common/test/debug
*
* initialize:
*
      do ii = 5,14
         cell(ii,ncell) = 0.0
      enddo
*
* calculate multipole moments for this cell
*
      do ii = i1,i2
         xdum  = lensdata(1,ii) - cell(1,ncell)              
         ydum  = lensdata(2,ii) - cell(2,ncell)
         xdum2 = xdum *xdum
         ydum2 = ydum *ydum
         xdum3 = xdum *xdum2
         ydum3 = ydum *ydum2
         xdum4 = xdum2*xdum2
         ydum4 = ydum2*ydum2
         xdum5 = xdum2*xdum3
         ydum5 = ydum2*ydum3
         xdum6 = xdum3*xdum3
         ydum6 = ydum3*ydum3
*
*    quadrupole terms:
*
         cell(5, ncell) = cell(5, ncell)
     &      + lensdata(3,ii)*( xdum2 - ydum2 )
         cell(6, ncell) = cell(6, ncell)
     &      + lensdata(3,ii)*2 * xdum*ydum
*
*    octopole terms:
*
         cell(7, ncell) = cell(7, ncell)
     &      + lensdata(3,ii)*( xdum3 - 3*xdum*ydum2 )
         cell(8, ncell) = cell(8, ncell)
     &      + lensdata(3,ii)*( 3*xdum2*ydum - ydum3 )
*
*    16-pole terms:
*
         cell(9, ncell) = cell(9, ncell)
     &      + lensdata(3,ii)*( xdum4 - 6*xdum2*ydum2 + ydum4 )
         cell(10,ncell) = cell(10,ncell)
     &      + lensdata(3,ii)*4 * ( xdum3*ydum - xdum*ydum3 )
*
*    32-pole terms:
*
         cell(11,ncell) = cell(11,ncell)
     &      + lensdata(3,ii)*(xdum5-10*xdum3*ydum2+5*xdum*ydum4)
         cell(12,ncell) = cell(12,ncell)
     &      + lensdata(3,ii)*(5*xdum4*ydum-10*xdum2*ydum3+ydum5)
*
*    64-pole terms:
*
         cell(13,ncell) = cell(13,ncell)
     &      + lensdata(3,ii)*(xdum6-15*xdum4*ydum2+15*xdum2*ydum4-ydum6)
         cell(14,ncell) = cell(14,ncell)
     &      + lensdata(3,ii)*(6*xdum5*ydum-20*xdum3*ydum3+6*xdum*ydum5) 
      enddo
*
      return
      end
