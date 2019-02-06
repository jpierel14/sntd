**********************************************************************  
**********************************************************************  
      subroutine input(debug,nray,arand,sigmas,sigmac,gamma,eps,minmass
     &,maxmass,power,pixmax0,pixminx,pixminy,pixdif,fracpixd
     &,ipix,ipix1,pix,pix1,times,eps2,powerp1
     &,powerp2,powerp3,aaa1,aaa2,string1,sigma,arand0,cr1,cr2,agam,bgam
     &,boa,bmsoams,ampth,pixmaxx,pixmaxy,facipix,pixdifh,rstars
     &,avcell,avlens,rayshot,lensmax,levelmax,pixlen,distcom2
     &,levellen,levelcel,nlens,factor1,factor2,massaver,raydist0
     &,raydist0sq,raydist1,raydist2,jxmin,jxmax,jymin,jymax
     &,stafac,del,raydifx,rayminx,rayminy,raymaxx,raymaxy,xur,yur,jobnu
     &,month,day,year,string)
*
* reads input data for microlensing program from file 'input',
* job number from file 'jobnu' and initializes some numbers
*
*
* input(parameters from main program):
*                ipix,ipix1,lensmax,levelmax,factor1,factor2
*
* output: 
*  	data read:
*		debug,nray,arand,sigmas,sigmac,gamma,eps,	
*               minmass,maxmass,power,pixmax0,pixminx,pixminy,
*               pixdif,fracpixd
*
*  	data initialized:
*                pix(),pix1(),times(20),
*                powerp1,powerp2,powerp3,aaa1,aaa2,string1,sigma,
*                arand0,cr1,cr2,eps2,agam,bgam,boa,bmsoams,ampth,
*                pixmaxx,pixmaxy,facipix,pixdifh,rstars,
*                nlens,avcell,rayshot,pixlen(2),distcom2,
*		 levellen,levelcel,massaver
*
* called in:
*           main
*
* use of subroutines:
*           none
*
* use of functions:
*           massav
*					October 31, 1991; J. Wambsganss
*
      implicit none
      external massav
      integer debug,nray,iwrite,factor2,factor1,ipix,ipix1,nlens
     &,idum,jdum,i1,jobnu,jxmax,jxmin,nrayx,nrayy,        lensmax
     &,levelmax,jymax,jymin,levellen(levelmax),levelcel(levelmax)
      integer*8 rayshot
*
      double precision arand,sigmas,sigmac,gamma,eps,minmass
     &,maxmass,power,pixmax0,pixminx,pixminy,fracpixd,pixdif
     &,powerp1,powerp2,powerp3,aaa1,aaa2,sigma,arand0,cr1,cr2
     &,eps2,agam,bgam,boa,bmsoams,ampth,pixdifh
     &,oopixfac,pixmaxx,pixmaxy,epsilon,pixdif0,facipix,denom
     &,rayminx,rayminy,rstars,raymaxx,raymaxy,raydist0sq,xur,yur
     &,raydist0,raydifx,raydist1,raydist2,del,distcom2(levelmax)
     &,xray,yray,massaver,massav,stamin,stamax,stadif,stafac
     &,avcell,avlens,dummy1,raydify,pixlen(2)
*
      integer*4 pix(ipix,ipix),pix1(ipix1,ipix1),month,day,year
      real times(20)
      character*3 string1, string2
      character*8 string
      common/leng/stamin,stadif
      data epsilon /1.0d-10/
*
* READING part:    
*
*
* determine job number (and add "1" to it):
*
      open(55,file='jobnum',status='old')
*
      read(55,'(a3)')string1
      read(string1,'(i3)')jobnu
*
* for increase of job number with every job:
*   comment out the next line, uncomment the
*	line below;
*
*
      write(string2,'(i3.3)')jobnu
c     write(string2,'(i3.3)')jobnu+ 1
      rewind 55 
      write(55,'(a3)')string2
      close(55)
*
*
*
*
      open(2,file='input',status='old')
*
      read(2,1000)arand,debug,sigmas,sigmac,gamma,eps,nray,minmass
     &,maxmass,power,pixmax0,pixminx,pixminy,pixdif,fracpixd,iwrite
      close(2)
*
* 
	 write(*,1094)jobnu,nray
	 write(*,1095)sigmas, sigmac, gamma
*
 1094 format(' jobnu = ',i3,' nray = ',i6)
 1095 format(' sigmas = ',f6.3,' sigmac = ',f6.3
     &	 ,' gamma = ',f6.3)
* 
*
*
*
* print out data just read (if iwrite > 0):
*
      if(iwrite.ne.0)then
         write(99,1000)arand,debug,sigmas,sigmac,gamma,eps,nray,minmass
     &	 ,maxmass,power,pixmax0,pixminx,pixminy,pixdif,fracpixd,iwrite
      endif
*
* INITIALIZING part:
*
      avcell   = 0.0
      avlens   = 0.0
      rayshot  = 0
*
      do idum = 1, ipix
        do jdum = 1, ipix
           pix(jdum,idum)  = 0
        enddo
      enddo
*
      do idum = 1, ipix1
         do jdum = 1, ipix1
            pix1(jdum,idum) = 0
         enddo
      enddo
*
      do i1 = 1,20
	times(i1) = 0.0e0
      enddo
*
      powerp1     = power + 1.0
      powerp2     = power + 2.0 
      powerp3     = power + 3.0 
*
      aaa2        = -(maxmass**powerp1)/powerp1
      aaa1        = -(minmass**powerp1)/powerp1 - aaa2
*
* total surface density:
*
      sigma = sigmas + sigmac
*
      if(abs(sigma-1.0).lt.epsilon)then
           write(44,*)' sigma = ',sigma,' !!! ---> stop 45 in INPUT' 
           stop 45
      endif
*
      arand0   = arand
      cr1      = 0.2
      cr2      = 0.3
      eps2     = eps*eps
      agam     = 1.0 + gamma
      bgam     = 1.0 - gamma
      boa      = abs(bgam/agam)
      bmsoams  = abs(  (bgam - sigma) / (agam - sigma)  )
      ampth    = 1.0d0/((1.0d0-sigma)**2-gamma**2)
*
* values for determination of stars:
*
      pixdif0  = 2.0*pixmax0
      oopixfac = float(ipix)/pixdif
      pixmaxx  = pixminx + pixdif
      pixmaxy  = pixminy + pixdif
      facipix  = float(ipix1)*oopixfac/(float(2*ipix))
      pixdifh  = 0.5*(pixmaxy - pixminy)
*
* lengths of shooting field in x-direction  
*    (proportionally: star fieldsdadd) 
*
      denom         = (1.0-sigma+gamma)
      if(abs(denom).gt.epsilon)then
         raymaxx    = (pixmaxx+fracpixd*pixdif0)/denom
         rayminx    = (pixminx-fracpixd*pixdif0)/denom
         if(rayminx.gt.raymaxx)then
	    dummy1  = rayminx
	    rayminx = raymaxx
	    raymaxx = dummy1    
	 endif
	 raydifx    = raymaxx - rayminx
         xray       = (pixmax0+5.0)/denom                       
      else
         write(44,*)' 1.0 - sigma + gamma = 0 ----> STOP 40 in INPUT !' 
         stop 40
      endif
*
* lengths of shooting field in y-direction  
*   (proportionally: star fieldsdadd) 
*
      denom         = (1.0-sigma-gamma)
      if(abs(denom).gt.epsilon)then
         raymaxy    = (pixmaxy+fracpixd*pixdif0)/denom
         rayminy    = (pixminy-fracpixd*pixdif0)/denom
         if(rayminy.gt.raymaxy)then
	    dummy1  = rayminy
	    rayminy = raymaxy
	    raymaxy = dummy1    
	 endif
	 raydify    = raymaxy - rayminy
         yray       = (pixmax0+5.0)/denom                       
      else
         write(44,*)' 1.0 - sigma - gamma = 0 ---> STOP 41 in INPUT !!!'
         stop 41
      endif
*
      rstars         = sqrt(xray**2 + yray**2) + 2                      
*
  222 continue
*
      massaver = massav()
*
*
*
      nlens    = max(nint(sigmas*rstars**2/massaver),1) 
*
      if(sigmas.lt.0.00001)then
         nlens = 1
      elseif(sigmas.lt.0.0015)then
         nlens = 2
      endif 
*
*
*
      stamin   =-rstars
      stamax   = rstars
      stadif   = stamax-stamin
      stafac   = stadif/float(ipix1)
*
      jymax    = nint(float(nray)*raymaxy/raydify)
      jymin    = nint(float(nray)*rayminy/raydify)
      jxmax    = nint(float(nray)*raymaxx/raydify)
      jxmin    = nint(float(nray)*rayminx/raydify)
      nrayx    = jxmax - jxmin
      nrayy    = jymax - jymin
*
* initialize vector with lengths of cells:
*        and vector with number of cells in certain levels:             
*
      do i1 = 1,levelmax
          distcom2(i1) = ( stadif/(2**i1) )**2
          levelcel(i1) = 0
          levellen(i1) = 0
      enddo
*
      if(nrayx.le.0)then
	 nrayx = 1
      endif
*
      raydist0   = raydifx/float(nrayx)
      raydist0   = raydify/float(nrayy)
      raydist1   = raydist0/float(factor1-1)
      raydist2   = raydist1/float(factor2)
      del        = 0.5*raydist1
      raydist0sq = raydist0**2
*
      pixlen(1)  = pixdif/float(ipix)
      pixlen(2)  = raydifx/float(ipix1)
*
* coordinates of  upper right corner of shooting rectangle:             
*
       xur  = (float(jxmax)+0.5)*raydist0 - 0.5*raydist2
       yur  = (float(jymax)+0.5)*raydist0 - 0.5*raydist2
*
       if(rstars**2.lt.xur**2+yur**2)then
    	  rstars  = sqrt(xur**2+yur**2) + epsilon 
 	  goto 222
       endif
*
      if(nlens.gt.lensmax)then
         write(44,*)' nlens = ',nlens,'   lensmax = ',lensmax          
         write(44,*)' lens number too large !!!!!'                     
         write(44,*)' STOP 1    in INPUT'
         stop 1
      endif
*
 1000 format(f10.3,/i10,4(/f10.3),/i10,2(/f10.6),6(/f10.3),/i10)
*
	return
	end
