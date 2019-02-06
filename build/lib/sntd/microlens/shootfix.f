**********************************************************************
**********************************************************************
      subroutine shootfix(cellmax,lensmax,xxx,yyy,alphax,alphay,cell
     &  ,lensdata,clno,lensused,totused)
*
*  is called for rays with fixed cells  and WITHOUT  lenses !!!
*     determination of ALPHAX, ALPHAY, not XSI, ETA !!!
*
*
* input:
*        xxx,yyy  - coordinates of ray in lens plane (before shooting)
*        cell     - array with (c.m.-) coord., mass, moments of cells
*        cellmax  - maximal number of cells
*        lensdata - array with coord., mass, moments of single lenses
*        lensmax  - maximal number of lenses
*        clno     - vector with indices of cells (positive) and lenses
*                   (negative) to be used
*        totused  - total number of used lenses and cells
* output:
*        xsi,eta - coordinates of ray in source plane (after shooting)
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
      integer cellmax,lensmax,iii,indum,lensused,totused,debug
     & ,clno(lensmax)
      double precision xxx,yyy,cell(14,cellmax),modist2,alphax,alphay
     &,dx,dx2,dx3,dx4,dx5,dx6,dx7,dy,dy2,dy3,dy4,dy5,dy6,dy7,p1,p2,p3
     &,p4,p5,p6,p7,p8,p9,p10,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,oodist2
     &,oodist6,oodist8,oodist10,oodist12,oodist14,lensdata(3,lensmax)
*
      common/test/debug
*
      alphax = 0.0D0
      alphay = 0.0D0
*
*  main loop over all cells:
*
      do iii = lensused + 1,totused
*
           indum = clno(iii)
*
* use this cell:
*
                dx        = (xxx - cell(1,indum))
                dy        = (yyy - cell(2,indum))
                dx2       = dx*dx
                dy2       = dy*dy
                oodist2 =          1.0      / (dx2 +  dy2)
                modist2 =     cell(3,indum)*oodist2
                alphax    =  alphax + dx*modist2
                alphay    =  alphay + dy*modist2
*
*   calculate and add quadrupole terms:
*
                          p1       = cell(5,indum)
                          p2       = cell(6,indum)
                          dx3      = dx*dx2
                          dy3      = dy*dy2
                          oodist6  = oodist2*oodist2*oodist2
                          q1       = dx3          - 3*dx*dy2
                          q2       =     3*dx2*dy          - dy3
*
                          alphax   = alphax + oodist6 * (p1*q1 + p2*q2) 
                          alphay   = alphay + oodist6 * (p1*q2 - p2*q1) 
*
*   calculate and add octopole terms:
*
                          p3       = cell(7,indum)
                          p4       = cell(8,indum)
                          dx4      = dx2*dx2
                          dy4      = dy2*dy2
                          oodist8  = oodist2*oodist6
                          q3       = +dx4         -6*dx2*dy2        +dy4
                          q4       =     +4*dx3*dy          -4*dx*dy3   
*
                          alphax   = alphax + oodist8 * (p3*q3 + p4*q4) 
                          alphay   = alphay + oodist8 * (p3*q4 - p4*q3) 
*
*   calculate and add 16-pole terms:
*
                          p5       = cell(9,indum)
                          p6       = cell(10,indum)
                          dx5      = dx2*dx3
                          dy5      = dy2*dy3
                          oodist10 = oodist2*oodist8
           q5       = +dx5         -10*dx3*dy2           +5*dx*dy4
           q6       =     +5*dx4*dy           -10*dx2*dy3         +dy5
*
                          alphax   = alphax + oodist10 * (p5*q5 + p6*q6)
                          alphay   = alphay + oodist10 * (p5*q6 - p6*q5)
*
*   calculate and add 32-pole terms:
*
                          p7       = cell(11,indum)
                          p8       = cell(12,indum)
                          dx6      = dx3*dx3
                          dy6      = dy3*dy3
                          oodist12 = oodist2*oodist10
           q7       =+dx6   -15*dx4*dy2           +15*dx2*dy4     -dy6  
           q8       =+6*dx5*dy         -20*dx3*dy3           +6*dx*dy5  
*
                          alphax   = alphax + oodist12 * (p7*q7 + p8*q8)
                          alphay   = alphay + oodist12 * (p7*q8 - p8*q7)
*
*   calculate and add 64-pole terms:
*
                          p9       = cell(13,indum)
                          p10      = cell(14,indum)
                          dx7      = dx3*dx4
                          dy7      = dy3*dy4
                          oodist14 = oodist2*oodist12
           q9       = dx7   -21*dx5*dy2           +35*dx3*dy4 -7*dx*dy6
           q10      = 7*dx6*dy   -35*dx4*dy3   +21*dx2*dy5  - dy7
*
                        alphax   = alphax + oodist14 * (p9*q9 + p10*q10)
                        alphay   = alphay + oodist14 * (p9*q10- p10*q9) 
*
      enddo
*
      return
      end
