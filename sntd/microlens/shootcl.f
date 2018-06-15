**********************************************************************  
**********************************************************************  
      subroutine shootcl(cellmax,lensmax,stack,xxx,yyy,cell
     &  ,eps2,cellused,icell,indcell,ncell,nlens,indlens,lensdata
     &  ,lensused,clno,lno,raydist0sq,totused)
*
*  is called for every single ray, works with stack (to simulate
*  recursive calls)
*
*  no actual calculation of xsi/eta or alphax/alphay, here only 
*  the cell/lens structure is determined !!! (--> vector CLNO)
*
* indcell and ncell are necessary only for function pos !
* 	istack   - height of stack
*
*
* input:
*        stack    - vector with indices of cells to work with
*        xxx,yyy  - coordinates of ray in lens plane (before shooting)  
*        cell     - array with (c.m.-) coord., mass, moments of cells   
*        eps2     - square of accuracy paramter
*        cellused - number of cells used for this particular ray
*        lensused - number of lenses used for this particular ray
*        totused  - sum of cellused and lensused
*        icell    - vector with actual numbers of cells
*        cellmax  - maximal number of cells
*        i1       -
*        indcell  -
*        ncell  - total number of cells
*        imult    - parameter that decides if higher multipole moments  
*                   should be used (1-no, 2-quad,... 6-64pole)
*        raydist0sq - distance**2 of rays in level 1 (cells closer to 
*			ray than that should be resolved)
*
* output:
*        xsi,eta - coordinates of ray in source plane (after shooting)  
*
*
* called in:
*           main
*
* use of subroutines:
*           none
*
* use of functions:
*           pos
*					October 31, 1991; J. Wambsganss
*
      implicit none
      external pos
*
      integer cellmax,lensmax,idum,i1,j1,indum,icell(cellmax),pos,ncell
     &,cellused,istack,stack(cellmax),guessc,jdum,nlens,lensused,debug
     &,clno(lensmax),lno(lensmax),indlens(lensmax),indcell(cellmax)
     &,totused,j2
*
      double precision dum2,xxx,yyy,cell(14,cellmax),eps2,dx,dx2,dy,dy2
     &,lensdata(3,lensmax),raydist0sq,oodist2
*
      common/test/debug
*
*     initialization: cellused, stack, istack, clno
*
      cellused = 0
      lensused = 0
      istack   = 0
      guessc   = 2
*
      do i1 = 1,lensmax
         clno(i1) = 0
      enddo
*
      do i1 = 4,1,-1
        jdum =pos(i1,indcell,cellmax,ncell,guessc,indlens,lensmax,nlens)
        if(jdum.ne.0)then
          istack   = istack + 1
          stack(istack) = jdum
        endif
      enddo
*
  100 continue
      if(istack.le.0) then
*
* shooting work done:   stack is empty !  
*     separate lenses from cells, then return to main:
*
*
* 	put cells/lenses used for this test ray in different vectors:
*
	    totused = lensused + cellused
            j1 = 0                                             
            do i1 = 1,totused                              
               if(clno(i1).lt.0)then                         
                  j1      = j1 + 1                              
                  lno(j1) = -clno(i1)                      
               endif                                         
	    enddo
*
            j2 = j1                                            
            do i1 = 1,totused                              
               if(clno(i1).gt.0)then 
                  j2      = j2 + 1
                  lno(j2) =  clno(i1)
               endif                                         
	    enddo
*
         return
      endif
*
*POPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPO       
*******************> pop         (start)    (take last element of stack)
*POPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPO       
*
      indum  = stack(istack)
      istack = istack - 1
*
*POPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPO       
*******************> pop         (end)
*POPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPOPO       
*
      if(indum.eq.0)then
*
* this cell/lens does not exist    ?!?
*
         write(44,*)' in shoot:  indum <= 0 !!! should not happen !!'  
         write(44,*)' ----> STOP  12  in SHOOT'
         stop 12
      elseif(indum.lt.0)then
*
* exactly one lens in cell:
*
         indum     = -indum
*
         lensused = lensused + 1
         clno(lensused + cellused) = - indum
*
*   ------------------> pop
*
           goto 100
      else
*
* more than one lens:   check, if cell size is smaller than given value:
*
         dx        = (xxx - cell(1,indum))
         dy        = (yyy - cell(2,indum))
         dx2       = dx*dx
         dy2       = dy*dy
         oodist2   = 1.0 / (dx2 +  dy2)
         dum2      = cell(4,indum)*oodist2
*
         if(dum2.ge.eps2.or.raydist0sq*oodist2.gt.1.0d0)then
*
*  more than one lens in cell  and  cell is too large to consider 
*
* or
*
* distance to (center of mass of) cell is smaller than distance between
*   two rays in level 1:
*
*   divide in four subcells
*
            idum = indcell(indum)*4
            do i1 = 4,1,-1
               j1 = idum + i1
*
****************************************************> push    (start)   
*
               jdum = pos(j1,indcell,cellmax,ncell,guessc,indlens
     &		       ,lensmax,nlens)
               if(jdum.ne.0)then
                  istack = istack + 1
                  stack(istack) =  jdum
                  guessc = jdum - 1
               endif
*
*****************************************************> push    (end)    
*
	    enddo
            goto 100
         else
*
* use this cell if (linear-length-of-cell/distance-of-ray-to-c.m.)<eps:
*
             cellused = cellused + 1
             clno(lensused + cellused) =  indum
*
* ------------------> pop
*
             goto 100
          endif
      endif
*
      end
