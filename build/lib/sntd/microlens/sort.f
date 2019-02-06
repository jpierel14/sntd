**********************************************************************  
**********************************************************************  
      subroutine sort(n,ra,rb,narray)
*
* sorts lenses/cells according to the INDEX of each lens/cell,
* which contains information about the cell it is in in each 
* higher level
*
* (routine from NUMERICAL RECIPES, slightliy changed by JKW)
*
*
* input :
*    n = nlens/ncell         ---> size of array to be sorted
*    ra = indlens/indcell    ---> array to be sorted in increasing order
*    rb = lensdata/cell      ---> data to be rearranged in this order
*    narray = 3/14		 ---> "depth" of array rb
*
* output: 
*    ra = indlens/indcell    ---> array,  sorted in increasing order
*    rb = lensdata/cell      ---> data,  rearranged in this order
*
* called in:
*           numstar
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
      integer n,i,j,k,ir,narray,iarray,ra(n),rra
      double precision rb(narray,n),rrb(20)
*
*
      k=n/2+1
      ir=n
10    continue
        if(k.gt.1)then
          k=k-1
          rra=ra(k)
	  do iarray = 1,narray
             rrb(iarray)=rb(iarray,k)
	  enddo 
        else
          rra=ra(ir)
	  do iarray = 1,narray
             rrb(iarray)=rb(iarray,ir)
	  enddo
*
          ra(ir)=ra(1)
	  do iarray = 1,narray
             rb(iarray,ir)=rb(iarray,1)
	  enddo
*
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
	    do iarray = 1,narray
               rb(iarray,1)=rrb(iarray)
	    enddo
            return
          endif
        endif
        i=k
        j=k+k
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
	    do iarray = 1,narray
               rb(iarray,i)=rb(iarray,j)
	    enddo
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        ra(i)=rra
	do iarray = 1,narray
           rb(iarray,i)=rrb(iarray)
	enddo
      go to 10
      end
