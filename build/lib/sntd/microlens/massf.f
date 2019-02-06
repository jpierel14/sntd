**********************************************************************  
**********************************************************************  
      double precision function massf(dum1)
*
* sets masses according to certain massfunction: power, minmass, maxmass
*           uses random numbers;
*
*
*  INPUT:
*      dum1        = random number in interval {0,1}
*
*  OUTPUT:
*      massf(function value)   = weighted mass according to power law
*
*
*
* called in:
*           setstar
*
* use of subroutines:
*           none
*
* use of functions:
*           none
*					October 31, 1991; J. Wambsganss
*
      implicit none
      double precision dum1,dum2,minmass,maxmass,power,powerp1
     &	,powerp2,powerp3,aaa1,aaa2
*
      common/m/minmass,maxmass,power,powerp1,powerp2,powerp3,aaa1,aaa2
*
* maps random number in interval {0,1}  to interval {0.740,16.583}:
*
      dum2 = -powerp1 * (aaa1*dum1 + aaa2) 
      massf = dum2**(1.0/powerp1)          
*
* {0.740,16.583}  --->  { minmass, maxmass }
*
      if(massf.lt.minmass-1.0e-10)then
	 write(44,*)' MASSF:  too low ! massf, minmass:',massf,minmass
	 write(44,*)'  ----------> STOP 345 in  MASSF'
	 stop 345 
      elseif(massf.gt.maxmass+1.0e-10)then
	 write(44,*)' MASSF:  mass too high !! ----> STOP 346 in  MASSF'
	 stop 346 
      else 
	 return
      endif
*
      end
************************************************************************
************************************************************************
      double precision function massav()   
*
* calculates average mass for given massfunction, minmass, maxmass;
*
*
*  input:
*      (just numbers from common block /m/ )
*
*  output:
*      massav(function value)   = average mass according to power law
*					and upper and lower cut-off
*
* called in:
*           input
*
* use of subroutines:
*           none
*
* use of functions:
*           none
*					October 31, 1991; J. Wambsganss
*
      implicit none
      double precision dumfm,dumf,minmass,maxmass,power,powerp1,powerp2
     &	,powerp3,aaa1,aaa2
*
      common/m/minmass,maxmass,power,powerp1,powerp2,powerp3,aaa1,aaa2
*
      if(maxmass-minmass.lt.1.0e-10)then
	 massav   = minmass
      else
	 dumfm    =  ( maxmass**powerp2 - minmass**powerp2 ) / powerp2 
	 dumf     =  ( maxmass**powerp1 - minmass**powerp1 ) / powerp1
	 massav   = dumfm/dumf 
      endif
*
      if(massav.lt.minmass)then
	 write(44,*)' average mass too low -------> STOP 365 in  MASSF'
	 stop 365
      elseif(massav.gt.maxmass)then
	 write(44,*)' average mass too high -------> STOP 366 in  MASSF'
	 stop 366
      else 
	 return
      endif
*
      end
