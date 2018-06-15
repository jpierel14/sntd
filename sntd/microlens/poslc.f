**********************************************************************  
**********************************************************************  
      integer function poslc(ndum,indlc,lcmax,index)                    
*
* used to find actual number of lens OR cell                 
* with index 'ndum'  ( searches with bisection in ordered (!) vector)   
*
* input:                                                                
*      ndum      - index of cell whose actual number is wanted          
*      indlens/indcell   - vector with actual numbers of lenses/cells   
*      lensmax/cellmax   - maximal number of lenses/cells               
*      index     - total number of (ordered) lenses /cells up to now    
*
* output:                                                               
*      position  -  function value (positive: exists at this position;  
*                  negative: does note exist, but belongs at - position;
*                  zero: not possible)                                  
*
* called in:                                                            
*           pos
*
* use of subroutines:                                                   
*           none                                                        
*
* use of functions:                                                     
*           none                                                        
*					October 31, 1991; J. Wambsganss
*
      implicit none                                                     
      integer ndum,lcmax,index,lower,upper,middle,debug,indlc(lcmax)
      common/test/debug                                                
*
* bisection with bracketing values lower, upper:                        
*
      lower = 0                                                         
      upper = index + 1                                                 
   10 if(upper-lower.gt.1)then                                          
         middle = (lower + upper)/2                                   
         if(ndum.gt.indlc(middle))then                                
            lower = middle                                          
         else                                                         
            upper = middle                                          
         endif                                                        
         goto 10                                                      
      endif                                                             
*
      if(upper.le.index.and.ndum.eq.indlc(upper))then                   
         poslc  = upper                                               
      elseif(lower.gt.0.and.ndum.eq.indlc(lower))then                   
         poslc  = lower                                               
      else                                                              
         poslc  = - upper                                             
      endif                                                             
*
      return                                                            
      end                                                               
