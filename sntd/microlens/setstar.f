**********************************************************************  
**********************************************************************  
      subroutine setstar(ar,cr1,cr2,rstars,nlens,lensmax,lensdata
     &				,masstot,pixdifh)
*
* nlens stars are randomly placed within a circle with radius rstars
* individual masses are attached to each star; 
*
* special treatment (end of this routine) if kappa  = 0:  
*	       then nlens = 1, put near the origin; 
*              for nlens = 2: positions/masses fixed
*
*  input:
*      ar          = a seed number > 0       for subroutine rands
*      cr1, cr2    = auxillary seed numbers   for subroutine rands      
*      rstars      = radius of circle with stars
*      nlens       = number of lenses
*
*  output:
*
*      lensdata(1,i) = x-coordinate of lens #i
*      lensdata(2,i) = y-coordinate of lens #i
*      lensdata(3,i) = mass         of lens #i
*
* called in:
*           main
*
* use of subroutines:
*           rands
*
* use of functions:
*           massf
*					October 31, 1991; J. Wambsganss
*
      implicit none
      real         time1,time2,times(20)
      external rands, massf
      integer nlens,ilens,lensmax,debug,nlensin
*
      double precision lensdata(3,lensmax),masstot,massf,pixdifh,ar,cr1
     &,cr2,rstars,randx,randy,randm
*
      common/test/debug
      common/jkwtime/times
*
* calculate two random numbers uniform in the range 0 to 1
*
         time1 = 0.0d0
         time2 = 0.0d0
      nlensin   = 0
      masstot   = 0.0
*
      do ilens = 1,nlens
    1    continue
*
	 call cpu_time (time1)
         call rands(ar,cr1,cr2,randx)
         call rands(ar,cr1,cr2,randy)
         call rands(ar,cr1,cr2,randm)
	 call cpu_time (time2)
         times(11)  = times(11) + (time2-time1)
*
         randx=randx*2.0D0-1
         randy=randy*2.0D0-1
*
* if position is too far from the center then do not use it   
*                 ( not inside the given circle !)
*
        if((randx**2+randy**2).gt.1.0d0)go to 1
*
* now the pair (randx,randy) fills uniformly a circle with radius 1;
* expand the random numbers so that they cover a circle with radius  
* "rstars";  determine random mass inside mass interval, according
* to power law:
*
         lensdata(1,ilens) = rstars*randx
         lensdata(2,ilens) = rstars*randy
         lensdata(3,ilens) = massf(randm)
*
	 masstot = masstot + lensdata(3,ilens)
	 if(abs(lensdata(1,ilens)).le.pixdifh.and.
     &	    abs(lensdata(2,ilens)).le.pixdifh)then
      		nlensin     = nlensin + 1
  	 endif
      enddo
*
*
*
*
*
*
*
*
* if there is only one lens: set it at the origin of the coordinates:
*
*      if there are two lenses: fix there position symmetrically 
*                    mass  2:1 
*
*
*
      if(nlens.eq.1)then
         lensdata(1,nlens) = 1.0e-10
         lensdata(2,nlens) = 1.0e-10
         lensdata(3,nlens) = 1.0e0
      elseif(nlens.eq.2)then
         lensdata(1,    1) = -0.92000000e0
         lensdata(2,    1) =  0.0000100e0
         lensdata(3,    1) =  2.000000e0
         lensdata(1,nlens) =  0.92000000e0
         lensdata(2,nlens) = -0.0000100e0
         lensdata(3,nlens) =  0.200000e0
      endif
*
      return
      end
