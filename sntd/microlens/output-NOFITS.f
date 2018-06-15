*********************************************************************  
**********************************************************************  
      subroutine output(ampth,rayshot,avlens,avcell,levelmax,levelcel
     &,levellen,highle,nlens,rstars,ncell,xur,yur,lensdata,lensmax
     &,indlens,cell,cellmax,indcell,pixlen,raydist0,raydist1,raydist2
     &,nray,nrayx,nrayy,rayminx,raymaxx,rayminy,raymaxy,jxmin,jxmax
     &,jymin,jymax,raydif,boa,bmsoams,factor1,factor2,jobnu,arand,debug
     &,sigmas,sigmac,gamma,eps,minmass, maxmass,power,massav,masstot
     &,pixmax0,pixminx,pixminy,pixdif,fracpixd,pix,pix1,ipix,ipix1,pixa
     &,pixb,month1,day1,year1,string1,stringno)
*
* does all the writing in a normal run
* 
*
* input :
*	(all data to be printed)
*
* output:
*	(none, except printing data)
*
* called in:                                                            
*           main 
*
* use of subroutines:                                                   
*           calcpix,calcpix2
*
* use of functions:                                                     
*           none                                                        
*					October 31, 1991; J. Wambsganss
*
      implicit none
      real times(20)
      integer nlens,ncell        ,cellmax,indcell(cellmax),pixhigh
     &,pixlow,ipix,debug
     &,                 rayth,levelmax,levelcel(levelmax),ipix1,i1,i2
     &,levellen(levelmax),highle,jxmin,jxmax,jymin,jymax,jobnu,nrayx,iii
     &,nray,nrayy,j1,j2,lensmax,indlens(lensmax),factor1,factor2
     &,maxpix1,ipix125,ipix175 
      integer*8 rayshot,raysarr ,rayslost
*
      integer*4 pix(ipix,ipix),pix1(ipix1,ipix1),month1,day1,year1
     &,month2,day2,year2
      integer*2 pixa(ipix,ipix),pixb(ipix1,ipix1)
*
      double precision rayperpi,rayamp1,ampth,average,average1
     &,avlens,rstars,avcell,xur,yur,pixlen(2),raydist0,raydist1
     &,raydist2,rayminx,raymaxx,rayminy,raymaxy,raydif,boa,bmsoams
     &,arand,sigmas,sigmac,gamma,eps,minmass
     &,maxmass,power,massav,masstot,pixmax0,pixminx,pixminy
     &,pixdif,fracpixd,cell(14,cellmax),lensdata(3,lensmax),raylev1
*
      character*8 string1,string2
      character*3 stringno
      character*2 string_vec(3) 
      integer*4 ivec1(3),ivec2(3)
      integer values(8)
*
      common/jkwtime/times
      integer ibin, ibin_max,idum
      parameter(ibin_max=400)
      integer mag_bin(0:ibin_max)
*
      integer iray_bin_max
      parameter(iray_bin_max = 100000)
      integer ray_bin(0:iray_bin_max)
*
* 
        integer   status, unit, blocksize, bitpix, naxis, naxes(2)
        integer   group, fpixel, nelements
               character filename*80      
               logical simple, extend
*


*
* determine day, month, year and system time (at end of job) : 
*
         call date_and_time(VALUES=values)
         day2   = values(3)
         month2 = values(2)
         year2  = values(1)
*
         write(string_vec(1),'(i2.2)')values(5)
         write(string_vec(2),'(i2.2)')values(6)
         write(string_vec(3),'(i2.2)')values(7)
         string2 = string_vec(1)//':'//string_vec(2)//':'//string_vec(3)
*
*
* determine pixels with  highest and lowest number of rays:
*
*
*
      do j1 = 0,iray_bin_max
	 ray_bin(j1) = 0
      enddo
*
      pixhigh  = 0
      pixlow   = 999999999
      raysarr  = 0
*
      do j2 = 1,ipix
	 do j1 = 1,ipix
	    iii  = pix(j1,j2)
	    raysarr = raysarr + iii
 	    if(pixhigh.lt.iii)then
               pixhigh = iii          
 	    elseif(pixlow.gt.iii)then
               pixlow = iii          
 	    endif
*
	    iii  = min(iii,iray_bin_max)
	    ray_bin(iii) = ray_bin(iii)  + 1
cw	    write(*,*)' j1,j2,iii,ray_bin(iii): ',j1,j2,iii,ray_bin(iii)
*
       	 enddo
      enddo
*
*
* determine number of lost rays, average number of ray per pixel, ...
*
      rayslost = rayshot - raysarr 
      rayth  = (jxmax-jxmin+1)*(jymax-jymin+1)*(factor1-1)**2*factor2**2
      raylev1  = float( (jxmax-jxmin+1)*(jymax-jymin+1) )
      rayperpi = float(raysarr)/float(ipix*ipix)
      rayamp1  = (pixlen(1)/raydist2)**2 
*
      avlens   = avlens/raylev1
      avcell   = avcell/raylev1
*
* opening files to write on:
*
      open(2,file='dat.'//stringno,status='unknown')
*
      write(2,1010)jobnu,day1,month1,year1,string1,day2,month2,year2
     &      ,string2,times(8),times(8)/60.0,times(8)/3600.0 
      write(2,1020)nlens,ncell,100.0*float(ncell)/float(nlens)
     &  ,avlens,avcell,times(7)/rayshot,rayperpi,pixhigh
     &  ,pixlow,rayamp1,rayperpi/rayamp1
     &   ,ampth,rayshot,raysarr,rayslost
      write(2,1030)arand,debug,sigmas,sigmac,gamma,eps,nray
     &     ,minmass, maxmass,power,pixmax0,pixminx,pixminy,pixdif
     &	   ,fracpixd
      write(2,1035) ipix,factor1,factor2
      write(2,1040)boa,bmsoams,massav,masstot,raydif
     &            ,jxmin,jxmax,jymin,jymax                           
     &            ,rayminx,raymaxx,rayminy,raymaxy
     &            ,raydist0,raydist1,raydist2,(pixlen(j1),j1=1,2)
      write(2,1045)rstars,xur,yur,rayth,avlens+avcell
     &  ,(avcell+avlens)/float(nlens),highle,indlens(nlens)
      if(times(7).le.0.0)then
         times(7) = 1.0
      endif
      write(2,1050)times(1),times(11),(times(j1),j1=2,3)
     &	   ,(times(j1),j1=12,15),times(7)
     &     ,(times(j1),100.0*times(j1)/times(7),j1=4,6)
     &     ,times(8)
*
      write(2,1060)nlens,' lenses '
      if(nlens.ge.10)then
           do j1 = 1,5
              write(2,1070)j1,indlens(j1),(lensdata(j2,j1),j2=1,3)
           enddo
           write(2,*)
           do j1 = nlens-4,nlens
              write(2,1070)j1,indlens(j1),(lensdata(j2,j1),j2=1,3)
           enddo
      else
           do j1 = 1,nlens
              write(2,1070)j1,indlens(j1),(lensdata(j2,j1),j2=1,3)
           enddo
      endif
*
      write(2,1060)ncell,' cells '
      if(ncell.ge.10)then
           do j1 = 1,5
              write(2,1070)j1,indcell(j1),(cell(j2,j1),j2=1,3)
           enddo
           write(2,*)
           do j1 = ncell-4,ncell
              write(2,1070)j1,indcell(j1),(cell(j2,j1),j2=1,3)
           enddo
      else
           do j1 = 1,ncell
              write(2,1070)j1,indcell(j1),(cell(j2,j1),j2=1,3)
           enddo
      endif
      write(2,1080)(j1,levelcel(j1),levellen(j1),j1=1,highle)
*
* writing file with number of rays per pixel:
*
c     open(3,file='data/CONUF'//stringno,status='new'
c    &   ,form='unformatted')
c     write(3)pix,pix1
c     close(3)
*
* now doing the "scaled" and formatted output for the IRIS:
*
      ampth    = abs(ampth)
      average  = rayamp1
      average1 = (float(ipix)/float(ipix1))**2*average*2.0   **2
*
      call calcpix(pix,ipix,ampth,rayamp1)
      call calcpix1(pix1,ipix1,ampth,average1)
*
      do i2=1,ipix
         do i1=1,ipix
 	    pixa(i1,i2) = pix(i1,i2)
 	 enddo
      enddo
*
      do i2=1,ipix1
         do i1=1,ipix1
 	    pixb(i1,i2) = pix1(i1,i2)
 	 enddo
      enddo
*
*
      maxpix1 = -1
      do i2=1,ipix1
         do i1=1,ipix1
            pixb(i1,i2) = pix1(i1,i2)
            maxpix1 = max(maxpix1,pixb(i1,i2))
         enddo
      enddo
*
      ipix125 = 0.25*ipix1
      ipix175 = 0.75*ipix1
      do i2=ipix125+1,ipix175
         pixb(ipix125,i2) = maxpix1
         pixb(ipix175,i2) = maxpix1
         pixb(i2,ipix125) = maxpix1
         pixb(i2,ipix175) = maxpix1
      enddo
*
*
*
*
      open(12,file='IRIS'//stringno,status='unknown',form='unformatted')
      write(12)pixa,pixb
      close(12)
*
*
*
*120918: in Cargese: error corrected, changed!
*
      do i1=1,ibin_max
cerror      mag_bin(ibin)  = 0
            mag_bin(i1)  = 0
      enddo

*
*
*
      do i2=1,ipix
         do i1=1,ipix
c	    ibin  = 100+100*log10(float(pixa(i1,i2))/rayamp1)
 	    ibin  =   100+100.0*(0.4*float(pixa(i1,i2)-1024)/256.0)
 	    ibin  = min(ibin,ibin_max)
 	    mag_bin(ibin)  = mag_bin(ibin)  +1
 	 enddo
      enddo
*
      idum = 0
*
*
*
* print output in FORMATTED file, for transfer to SUNs:
*
c     open(12,file='ISIS'//stringno,status='unknown')
c         write(12,'((20i4))')pixa,pixb
c	  close(12)

*
*
*
*
c     open(56,file='ray_stat_'//stringno,status='unknown')
c     do j1 = 0,iray_bin_max
cc	 write(*,*)' j1,ray_bin(j1): ', j1,ray_bin(j1)
c	 if(ray_bin(j1).ne.0)then
c	    write(56,*)j1,ray_bin(j1)
c	 endif
c     enddo
*
*
*
      return
*
* different formats: 
*
 1010 format(i20   ,' jobnu:   number of job'
     &      / 20x,'date and time of start:' i4,'.',i2,'.',i4,'    ',a8
     &      / 20x,'date and time of stop: ' i4,'.',i2,'.',i4,'    ',a8
     &      /'CPU-time: ',f12.2,' sec = ',f10.2,' min = ',f8.2,' h')
 1020 format(/' important parameters:'
     &/,i20   ,' nlens:      total number of lenses within rstars'  
     &/,i20   ,' ncell:      total number of cells'
     &/,f20.2 ,'% ncell/nlens: ratio in percent'
     &/,f20.3 ,' avlens:     average number of lenses used per ray'
     &/,f20.3 ,' avcell:     average number of cells  used per ray'
     &/,f20.9 ,' CPU-time/ray: shooting time (CPU-sec) per ray' 
     &/,f20.3 ,' rayperpi:   average number of rays per pixel'
     &/,i20   ,' pixhigh:    highest number of rays per pixel'
     &/,i20   ,' pixlow:     lowest number of rays per pixel'
     &/,f20.3 ,' rayamp1:    number of rays for amplification 1'
     &/,f20.3 ,' ampav:      (numerical)   average amplification'
     &/,f20.3 ,' ampth:      (theoretical) average amplification'
     &/,i20   ,' rayshot:    total number of rays shot in level 3'  
     &/,i20   ,' raysarr:    number of rays arrived in square'
     &/,i20   ,' rayslost:   number of lost rays')
 1030 format(/' INPUT parameters:'
     &/,f20.3 ,' arand:   input for random number generator'
     &/,i20   ,' debug:   parameter for debugging the program       '
     &/,f20.3 ,' kappa:   surface mass density in compact objects   '
     &/,f20.3 ,' kappa_s: surface mass density in smooth matter     '
     &/,f20.3 ,' gamma:   global shear'
     &/,f20.3 ,' eps:     accuracy parameter, (0 <= eps <= 1)'
     &/,i20   ,' nray:    number of rays per row in level 1'
     &/,f20.3 ,' minmass: lower cutoff of mass spectrum'
     &/,f20.3 ,' maxmass: upper cutoff of mass spectrum'
     &/,f20.3 ,' power:   exponent of mass spectrum (Salpeter: 2.35)'
     &/,f20.3 ,' pixmax0: size of field for distribution of stars   '
     &/,f20.3 ,' pixminx: left border of receiving field'
     &/,f20.3 ,' pixminy: lower border of receiving field'
     &/,f20.3 ,' pixdif:  size of receiving field'
     &/,f20.3 ,' fracpixd:  fraction added ')
 1035 format(
     &/,i20   ,' ipix:    size of pixel matrix  IRISxxx'
     &/,i20   ,' factor1: multiplier: # rays from level 1 to level 2'
     &/,i20   ,' factor2: multiplier: # rays from level 2 to level 3')
 1040 format(/' OUTPUT parameters:'
     &/,f20.3 ,' boa:     b / a  :     (1-gamma) / (1+gamma)'
     &/,f20.3 ,' bmsoams: b-s/a-s: (1-gamma-sigma) / (1+gamma-sigma)'
     &/,f20.3 ,' massav:  average mass of lenses (in solar masses)'
     &/,f20.3 ,' masstot: total mass in all lenses (in solar masses)'
     &/,f20.3 ,' raydif:  length of shooting region (level 1)'
     &/,i20   ,' jxmin:   number of leftmost  ray in level 1'
     &/,i20   ,' jxmax:   number of rightmost ray in level 1'
     &/,i20   ,' jymin:   number of lowest    ray in level 1'
     &/,i20   ,' jymax:   number of highest   ray in level 1'
     &/,f20.3 ,' rayminx  coordinate of leftmost  ray in level 1'
     &/,f20.3 ,' raymaxx  coordinate of rightmost ray in level 1'
     &/,f20.3 ,' rayminy  coordinate of lowest    ray in level 1'
     &/,f20.3 ,' raymaxy  coordinate of highest   ray in level 1'
     &/,f20.6 ,' raydist0: distance of adjacent rays in level 1'
     &/,f20.6 ,' raydist1: distance of adjacent rays in level 2'
     &/,f20.6 ,' raydist2: distance of adjacent rays in level 3'
     &/,f20.6 ,' pixlen(1): pixel length in units of Einstein radii'
     &/,f20.6 ,' pixlen(2): pixel length in lens array PIXLENS')
 1045 format(
     &  f20.3 ,' rstars:    radius of circle that comprises stars'  
     &/,f20.6 ,' xur:       coordinates of leftmost ray in level 3'
     &/,f20.6 ,' yur:       coordinates of lowest   ray in level 3'
     &/,i20   ,' rayth:     number of rays expected to be shot'
     &/,f20.3 ,' avtot:     average # of cells+lenses  used per ray'
     &/,f20.6 ,' avtot/nlens: avtot divided by total # of lenses'
     &/,i20   ,' highle:    highest level used for cell hierarchy'
     &/,i20   ,' indlens(nlens) number of highest lens')
 1050 format(/' CPU time used for different parts of the program:'
     &/,f20.3 ,' times(1):   time in CPU-sec for subroutine setstar'
     &/,f30.3 ,' times(11): CPU-sec for subroutine rands'
     &/,f20.3 ,' times(2):   time in CPU-sec for subroutine numstar'
     &/,f20.3 ,' times(3):   time in CPU-sec for subroutine setcell'
     &/,f30.3 ,' times(12): CPU-sec for subroutine setmult'
     &/,f30.3 ,' times(13): CPU-sec for     setup    '
     &/,f30.3 ,' times(14): CPU-sec for sorting lenses'
     &/,f30.3 ,' times(15): CPU-sec for sorting cells'
     &/,f20.3 ,' times(7):   time in CPU-sec for    shooting'
     &/,f30.3 ,' times(4): CPU-sec for shootcl  (',f5.2,'%)'
     &/,f30.3 ,' times(5): CPU-sec for shootfix (',f5.2,'%)'
     &/,f30.3 ,' times(6): CPU-sec for shotayco (',f5.2,'%)'
     &/,f20.3 ,' times(8):   time in CPU-sec whole run') 
 1060 format(//,i10,a8'  after sorting: ',/
     & 9x,'number',15x,'index      x-coord       y-coord      mass')
 1070 format(i15,i20,2f14.6,f9.2)
 1080 format(//' distribution of cells and lenses in different levels:'
     &      ,/,13x,'level',14x,'cells',15x,'lenses',/(15x,i3,2i20))
 1090 format('  i3   raysarr  rayslost  pixhigh   pixlow  rayperpi',
     &  '     ampth ampactual')
 1100 format(i3,2i15,i8,i4,3f9.3)
      end                                                               
*
**********************************************************************  
**********************************************************************  
      subroutine calcpix(pix,ipix,amp,average)
*
* calulates "color" values of pix, pix1:
*
*
* mapping of pixel field on integers: 
*
*  average amplification  --->   1024  (irismed)
*   one   mag brighter than average  ---> 1024+1*256 = 1280
*   two   mag brighter than average  ---> 1024+2*256 = 1536
*   three mag brighter than average  ---> 1024+3*256 = 1792
*   one   mag dimmer   than average  ---> 1024-1*256 =  768
*   two   mag dimmer   than average  ---> 1024-2*256 =  512
*   three mag dimmer   than average  ---> 1024-3*256 =  256
*
*					October 31, 1991; J. Wambsganss
*
      implicit	none
      integer  pixmin,pixmax,ipix,i,j,irismed
      integer*4 pix(ipix,ipix)
      double precision factor,amp,average,zero,log10pix
*
* zero: logarithm of (number of rays corresponding to 
*       amplification 1  times   theoretical amplification)
*
      irismed = 1024
      factor  = 256.0
      zero    = log10(average*amp)
      do j = 1,ipix
         do i = 1,ipix
            if(pix(i,j).le.0)then
	       pix(i,j) = 1
	    endif
            log10pix = log10(float(pix(i,j)))
            pix(i,j) =irismed+2.5*(log10pix-zero)*factor
*
*020102:  produce STEPS in color distribution:
*
c           pix(i,j) =irismed+nint(2.5*(log10pix-zero))*factor
*020403:
c           if(pix(i,j).lt.256)then
c              pix(i,j) =256
c	    endif
         enddo
      enddo 
*
* determine min and max: 
*
         pixmax  = 0
         pixmin = 99999
         do j = 1,ipix
            do i = 1,ipix
 	       pixmin = min(pix(i,j),pixmin)
 	       pixmax = max(pix(i,j),pixmax)
 	    enddo
 	 enddo 
*
      return
      end
**********************************************************************  
**********************************************************************  
	subroutine calcpix1(pix,ipix,amp,average)
*
* calulates "color" values of pix, pix1:
*
*
* mapping of pixel field on integers: 
*
*  average amplification  --->   1024  (irismed)
*   one   mag brighter than average  ---> 1024+1*256 = 1280
*   two   mag brighter than average  ---> 1024+2*256 = 1536
*   three mag brighter than average  ---> 1024+3*256 = 1792
*   one   mag dimmer   than average  ---> 1024-1*256 =  768
*   two   mag dimmer   than average  ---> 1024-2*256 =  512
*   three mag dimmer   than average  ---> 1024-3*256 =  256
*
*					October 31, 1991; J. Wambsganss
*
  	implicit	none
	integer  pixmin,pixmax,ipix,i,j,irismed
	integer*4 pix(ipix,ipix)
	double precision factor, amp, average, zero,log10pix
*
* zero: logarithm of (number of rays corresponding to 
*       amplification 1  times   theoretical amplification)
*
        zero    = log10(average*amp)
	irismed = 1024
	factor  = 256.0
        do j = 1,ipix
           do i = 1,ipix
	     if(pix(i,j).le.0)then
		pix(i,j) = 1
	     endif
	     log10pix = log10(  float( pix(i,j) )  )
             pix(i,j) = irismed+2.5*(log10pix-zero)*factor
             if(pix(i,j).lt.256)then
	        pix(i,j)=256
	     endif
	   enddo
 	enddo 
*
* determine min and max: 
*
        pixmax  = 0
        pixmin = 99999
        do j = 1,ipix
           do i = 1,ipix
 	      pixmin = min(pix(i,j),pixmin)
 	      pixmax = max(pix(i,j),pixmax)
 	  enddo
 	enddo 
*
	return
	end
