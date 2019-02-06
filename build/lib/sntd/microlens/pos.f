**********************************************************************  
**********************************************************************  
      integer function pos(ndum,indcell,cellmax,ncell,guessc
     &                     ,indlens,lensmax,nlens)
*
* determines number of cell or lens in shoot  (new version: bisection)  
*     implemetation of "hunt" (--> num.rec.)
*
* input:
*      ndum      - index of cell whose actual number is wanted
*      indcell   - vector with actual numbers of cells
*      cellmax   - maximal number of cells
*      ncell   - total number of cells
*      guessc    - guess for actual cell number; accelerates search a lo
*                      (in case cell search is negative)
*      indlens   - vector with actual numbers of lenses
*      lensmax   - maximal number of lenses
*      nlens     - total number of lenses
*
* output(function value):
*        - positive: position of data for cell number "ndum" in vectors 
*          indcell, cell(i,.),...   (if cell with # ndum exists)        
*        - negative: position of data for lens number "ndum" in vectors 
*          indlens, lensdata(i,.),...   (if lens with # ndum exists)    
*        - zero     (if neither cell nor lens with # ndum does exist)   
*
* called in:
*           shootcl
*
* use of subroutines:
*           none
*
* use of functions:
*           poslc
*					October 31, 1991; J. Wambsganss
*
      implicit none
      external poslc
      integer cellmax,ncell,guessc,plens,poslc,lensmax,nlens,inc,lower
     &      ,upper,middle,ndum,indlens(lensmax),indcell(cellmax)
*
      if(guessc.lt.2)guessc = 2
*
* first check, if guess fits:
*
      if(ndum.eq.indcell(guessc))then
         pos = guessc
         return
*
* then check, if number of cell larger than largest value in indcell:   
*               (--> look for according lens)
*
      elseif(ndum.gt.indcell(ncell))then
         plens =  poslc(ndum,indlens,lensmax,nlens)
         if(plens.gt.0)then
            pos = -plens
         else
            pos = 0
         endif
         return
*
* then check, if number of cell less than smallest value in indcell:    
*               (--> return with zero)
*
      elseif(ndum.lt.indcell(1))then
         if(ndum.lt.indlens(1))then
            pos = 0
            return
         else
            plens =  poslc(ndum,indlens,lensmax,nlens)
            if(plens.gt.0)then
               pos = -plens
            else
               pos = 0
            endif
            return
         endif
      else
*
*  start with "HUNT"-procedure:   ( guess value: guessc )
*
         inc   = 1
         upper = guessc
         lower = guessc - 1
*
*    bracket searched value with  lower  and  upper:
*
    5    if(ndum.lt.indcell(lower))then
            upper  = lower
            lower = max(lower-inc,1)
            inc    = inc + inc
            goto 5
         elseif(ndum.gt.indcell(upper))then
            lower  = upper
            upper = min(upper+inc,ncell)
            inc    = inc + inc
            goto 5
         endif
*
* bisection with bracketing values lower, upper:
*
   10    if(lower.gt.upper)then
*
*  no cell found: look for lens with this (cell-) number:
*
            plens =  poslc(ndum,indlens,lensmax,nlens)
            if(plens.gt.0)then
               pos = -plens
            else
               pos = 0
            endif
            return
         else
*
*  next step in bracketing:
*
            middle = (lower + upper)/2
            if(ndum.eq.indcell(middle))then
               pos = middle
               return
            elseif(ndum.lt.indcell(middle))then
               upper = middle - 1
               goto 10
            else
               lower = middle + 1
               goto 10
            endif
         endif
      endif
      end
