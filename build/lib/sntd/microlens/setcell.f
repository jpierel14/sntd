**********************************************************************  
**********************************************************************  
      subroutine setcell(cellmax,lensmax,levelmax,cell,icell,indcell
     &  ,ncell,distcom2,lensdata,indlens,nlens,highle,levelcel,levellen)
*
* assigns new numbers to indlens(.) according to cell number in which 
* this lens is single by comparison of two neighbour lenses in 'old'
* array indlens(.)  (--> numstar.f) . Then sorts lenses according 
* to 'new' indlens(.). Determines number of lenses and cells per niveau
*
* input : lensdata(1:3,i)          coordinates and mass of lens #i
*         indlens(lensmax)         ordered (cell-)numbers of lenses
*         cellmax                  maximal number of cells
*         levelmax                 maximal number of levels
*         distcomp2(i)             size of cells in level  i
*         indcell(i)               actual number of cell with index i
*         ncell                    total number of cells
*         nlens                    total  number of lenses
*         ilens                    current number of lens
*
* output:
*         cell( 1:14,.),icell(.)   data of all cells
*	  highle		   level of highest lens
*	  levelcel(levelmax)	   number of cells per level
*	  levellen(levelmax)	   number of lenses per level
*				(cells with one lens in it)
*
* called in:
*           main
*
* use of subroutines:
*           setmult, sort
*
* use of functions:                                                     
*           niveau
*					October 31, 1991; J. Wambsganss
*
      implicit none
      external setmult,sort,niveau
c     real cputime,times(20)
      real         times(20)
      integer cellmax,lensmax,levelmax,icell(cellmax),ncell,nlens,ilens
     &,debug,j,lev,k,m,firstniv(30),jmax,fourpow(31),leveldif,index1
     &,index2,tempold(0:30),indlens(lensmax),i1,i2,j1,temp(0:30),highle
     &,indcell(cellmax),niveau,levelcel(levelmax),levellen(levelmax)
      double precision lensdata(3,lensmax),cell(14,cellmax),time1,time2
     &,distcom2(levelmax),times1,times2,tsetmult,fourpowi(30)
      common/test/debug
      common/jkwtime/times
*
*     initialize:
*
      jmax     = 0
      ncell    = 0
      leveldif = 0
      temp(0)  = 0
      tsetmult = 0.0
*
      do j = 1,levelmax
         firstniv(j) = 0
	 fourpow(j)  = 4**(j-1)
	 fourpowi(j) = 1.0/fourpow(j) 
      enddo
      fourpow(levelmax+1) = 0
*
      call cpu_time (time1)
*
      do 10 ilens = 1,nlens
*
*  loop over all lenses (compare lens #ilens with lens # ilens+1):
*
         index1 =  indlens(ilens)
	 index2 =  indlens(ilens+1)
*
	 do 20 j = 1,levelmax
*
*  loop over all levels (up to the first with differing numbers)
*
	    index1  =  index1  - i1*fourpow(levelmax+2-j)
	    index2  =  index2  - i1*fourpow(levelmax+2-j)
*
	    i1      =  index1  * fourpowi(levelmax+1-j)
	    i2      =  index2  * fourpowi(levelmax+1-j)
 	    temp(j) =  i1 
*
            if(i1.ne.i2)then
*
* the two lenses are in different cells in level j:
* 
               jmax     =  max(jmax,j)
	       if(j.lt.leveldif)then
*
*  level smaller than in case before: (-> create cells )
*
	          indlens(ilens) = 0
		  do k = 1,leveldif
		     indlens(ilens) = indlens(ilens)+
     &        		(tempold(k)+1) * fourpow(leveldif-k+1)
		  enddo
*
		  do k = 1,j-1         
		     tempold(k) = temp(k) 
	 	  enddo
		  tempold(j) = i2
*
*  create cells on level j and higher:
*
	          do lev = j,leveldif-1
		     ncell          = ncell + 1
		     if(ncell.gt.cellmax)then
                       write(44,*)' ncell =',ncell,' cellmax =',cellmax
                       write(44,*)' ncell > cellmax:  STOP 2 in SETCELL'
                       stop 2
                     endif 
		     indcell(ncell) =  0
 		     do m = 1,lev
	 	        indcell(ncell) = indcell(ncell) +
     &    		    (temp(m)+1) * fourpow(lev-m+1)
		     enddo
*
* setting up center of mass coordinates and total mass of this cell:
*
	             do m = firstniv(lev),ilens 
 			cell(1,ncell) = cell(1,ncell) 
     &				+ lensdata(1,m)*lensdata(3,m)
 			cell(2,ncell) = cell(2,ncell) 
     &				+ lensdata(2,m)*lensdata(3,m)
 			cell(3,ncell) = cell(3,ncell) 
     &				+ lensdata(3,m)
 		     enddo 
*
		     cell(1,ncell) = cell(1,ncell) / cell(3,ncell)
 		     cell(2,ncell) = cell(2,ncell) / cell(3,ncell)
 		     cell(4,ncell) = distcom2(lev)
		     icell(ncell) = ilens - firstniv(lev) + 1
*
                     call cpu_time (time1)
		     call setmult (lensdata,lensmax,nlens
     &				,cell,cellmax,ncell,firstniv(lev),ilens)
                     call cpu_time (time2)
		     tsetmult = tsetmult + times2-times1
		     times(12) = times(12) + times2-times1
*
		     firstniv(lev) = -1 
	          enddo 
	       else
*
* level bigger than in case before:
*
	          indlens(ilens) = 0

		  do k = 1,j
 		     indlens(ilens) = indlens(ilens)+
     &        	       	(temp(k)+1) * fourpow(j-k+1)
		  enddo
		  do k = 1,j-1
		     tempold(k) = temp(k) 
	 	  enddo
	          tempold(j) = i2
	       endif 
	       leveldif = j 
	       goto 10
	    else  
*
* 		the two lenses are in the SAME cell in level j:
* 
	       if(j.eq.levelmax)then
		  write(44,1050)ilens,indlens(ilens),
     &				     indlens(ilens+1),j,i1,i2
	 	  stop 154
	       endif
* 
	       if(firstniv(j).le.0)then
	          firstniv(j) = ilens
	       endif
* 
 	    endif
* 
   20    continue
   10 continue
*
	time1 = time2
        call cpu_time (time2)
        times(13) = time2-time1 - times(12) 
        time1 = time2
*
	if(nlens.gt.1)then
            call sort(nlens,indlens,lensdata,3)
	endif
*
        call cpu_time (time2)
        times(14) = times(14) + (time2-time1)
        time1 = time2
*
	if(ncell.gt.1)then
        	call sort(ncell,indcell,cell,14)
	endif
*
*     determine lens of highest level used:
*
      if(nlens.gt.0)then
           highle = niveau(indlens(nlens))
      else
           highle = 0
      endif
*
*     count  cells of different levels
*
      do i1= 1,ncell
           j1 = niveau(indcell(i1))
           levelcel(j1) = levelcel(j1) + 1
      enddo
*
*     count  lenses of different levels
*
      do i1= 1,nlens
           j1           = niveau(indlens(i1))
           levellen(j1) = levellen(j1) + 1
      enddo
*
        call cpu_time (time2)
      times(15) = times(15) + (time2-time1)
*
 1050 format(' in SETCELL: need more levels, levelmax too small !!! ',/,
     & ' ilens,indlens(ilens),indlens(ilens+1),j,i1,i2',/,6i20
     &	/'   ---------> STOP 154 in SETCELL')
*
      return
      end
