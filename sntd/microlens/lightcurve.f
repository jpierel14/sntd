	program readIRIS
*
*      reads   unformatted  (pixel)  data of 1000^2 pixels from file IRIS401
*       draws a line through two points, determines lightcurve
*	for circular source, radius sigma,  with flat profile 
*	   (Gaussian source with width sigma possible as well) 
*
*
*
*					Joachim Wambsganss, January 2008
*							jkw, 2/27/96
*
*
      implicit none
*
      integer ipix,ipix1,pixlow,i1,i2,iii,ix1,ix2,iy1,iy2,i0
     &	,ixxx,iyyy,pixmax,ilines,iamp,sigma,Rn
     &	,ipoints,ix,iy
      parameter(ipix=1000,ipix1=500)
      integer*2 pix(ipix,ipix),pix1(ipix1,ipix1)
      double precision aaa,x1,x2,y1,y2,slope, alpha,sinalpha,cosalpha
     &	,x0,y0,x_end,y_end,x_start,y_start,x_diff,y_diff,xxx,yyy,value
     &  ,a,b
      real  pix_real(ipix,ipix)
      character*3 string3
      character*1 string1
      character*80 fileread, filewrite
      parameter(ipoints = 505)
      real  x(ipoints),y(ipoints), xx(ipoints),yy(ipoints)

*       My meddling
	open(1,file='source'status='old')
	read(1,1000)startLocx,startLocy,endLocx,endLocy,radius
	
*       
* name of magnification pattern to be read:
*
	fileread  = 'IRIS001'
	
*
* name of magnification pattern to be written on:
*
      filewrite = 'IRIS001-track'
*
*     
* position of start point and end point of lightcurve:
*
*
	ix1 = startLocx
	iy1 = startLocy	
	ix2 = endLocx
	iy2 = endLocy
*
* source size in pixels (Gaussian width)
*
       sigma = radius
       Rn=20       	
       a=-5
       b=0	
*
         x1 = dble(ix1)
         x2 = dble(ix2)
         y1 = dble(iy1)
         y2 = dble(iy2)
*
*
*
*
      write(*,*)' reading from file ',fileread
      open(3,file=fileread ,status='old',form='unformatted')
      read(3)pix,pix1
      close(3)
      write(*,*)'                           ...  done'
*
* determine maximum of magnification pattern:
*
       pixmax = 0
       do i1 = 1,ipix
          do i2 = 1,ipix
	      pixmax = max(pixmax,pix(i2,i1))
	     pix_real(i2,i1) = 10**(0.4*(float(pix(i2,i1)-1024)/256.0))
          enddo
       enddo
       write(*,*)' pixmax = ',pixmax
*
*
* determine points of line:
*
         slope = (y2-y1)/(x2-x1)
         write(*,*)' slope = ',slope
*
         alpha = atan(slope)
         write(*,*)' alpha = ',alpha
         sinalpha = sin(alpha)
         cosalpha = cos(alpha)
*
* points for x = 1    and    y = 1:
*       
	 if(slope.ne.0) then
	    y0  =   y1 + (1.0d0-x1)*slope
	    x0  =   x1 + (1.0d0-y1)/slope
	    write(*,*)' y0 = ',y0
	    write(*,*)' x0 = ',x0
	 else
	    y0=y1
	    x0=x1
*
*
*
         if(x0.ge.1.0d0.and.x0.le.dble(ipix))then
            x_start = x0
            y_start = 1.0d0
         elseif(y0.ge.1.0d0.and.y0.le.dble(ipix))then
            x_start = 1.0d0
            y_start = y0
         endif

         write(*,*)' x_start = ',x_start
         write(*,*)' y_start = ',y_start
*
*
*
* points for x = ipix    and    y = ipix:
*       
	 if(slope.ne.0) then
	    y0  =   y1 + (dble(ipix)-x1)*slope
	    x0  =   x1 + (dble(ipix)-y1)/slope
	    write(*,*)' y0 = ',y0
	    write(*,*)' x0 = ',x0
	 else
	    x0=x1
	    y0=y1
*       
*
*
         if(x0.ge.1.0d0.and.x0.le.dble(ipix))then
            x_end = x0
            y_end = y1 + (x_end-x1)*slope
         elseif(y0.ge.1.0d0.and.y0.le.dble(ipix))then
            y_end = y0
            x_end = x1 + (y_end-y1)/slope
         endif
         write(*,*)' x_end = ',x_end
         write(*,*)' y_end = ',y_end
     
         x_diff = x_end - x_start
         y_diff = y_end - y_start
         write(*,*)' x_diff = ',x_diff
         write(*,*)' y_diff = ',y_diff
*
*
*
* mark pixels along line:
*
* AND   determine light curve for Gaussian source:
*
*
*
         open(23,file='out_line',status='unknown')
c        do i0 = 0,2*ipix
         do i0 = -2*ipix,2*ipix
	    
	    xxx = x_start + i0*cosalpha
	    yyy = y_start + i0*sinalpha
	    ixxx = nint(xxx)
	    iyyy = nint(yyy)
	    if(ixxx.ge.1+3*sigma.and.ixxx.le.ipix-3*sigma.and.
     &	          iyyy.ge.1+3*sigma.and.iyyy.le.ipix-3*sigma)then
*
*
                 call gauss(sigma,Rn,a,b,xxx,yyy,value,pix_real,ipix)
*
*
*
ccc	         value = 10**(float(pix(ixxx,iyyy)-1024)/256.0)
	         write(23,1010)i0,xxx,yyy,pix(ixxx,iyyy),value,ixxx,iyyy
            endif
	 enddo
         do i0 = -2*ipix,2*ipix
	    xxx = x_start + i0*cosalpha
	    yyy = y_start + i0*sinalpha
	    ixxx = nint(xxx)
	    iyyy = nint(yyy)
	    if(ixxx.ge.1+3*sigma.and.ixxx.le.ipix-3*sigma.and.
     &	          iyyy.ge.1+3*sigma.and.iyyy.le.ipix-3*sigma)then
                  pix(ixxx,iyyy) = 0
                  pix(ixxx,iyyy) = 2000
                  pix(ixxx,iyyy) = 5000
            endif
	 enddo
 1010 format(i5,2f15.5,i8,f12.5,2i5)

*
*
*030711: indicate source in lower left corner:
*
*
         do i1 = -2*sigma, 2*sigma
            do i2 = -2*sigma, 2*sigma
		 if(i1*i1+i2*i2.le.2*sigma*2*sigma)then
	             pix(10+2*sigma+i1,10+2*sigma+i2) = 2*pixmax/3
		 endif
	    enddo
	 enddo

         do i1 = -sigma, sigma
            do i2 = -sigma, sigma
		 if(i1*i1+i2*i2.le.sigma*sigma)then
	             pix(10+2*sigma+i1,10+2*sigma+i2) = pixmax
		 endif
	    enddo
	 enddo



*
*
* convert x and y into pixels: 
* 		and add pixel values to magpat: 
*
      do i1 = 1,ipoints
        ix = ((x(i1) +1.1)*1000/2.2) + 1
        iy = ((y(i1) +1.1)*1000/2.2) + 1
c	write(*,*)' i1,x(i1),ix,y(i1),iy: ', i1,x(i1),ix,y(i1),iy
	if(ix.ge.1.and.ix.le.ipix.and.iy.ge.1.and.iy.le.ipix)then
c	   pix(ix,iy) = 5000
c	   pix(ix,1000-iy) = 5000
	   pix(ix,iy)      = pix(ix,iy)      + 500
	   pix(ix,1000-iy) = pix(ix,1000-iy) + 500
	else
c	   write(*,*)' i1,x(i1),ix,y(i1),iy: ', i1,x(i1),ix,y(i1),iy
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
*
      write(*,*)' filewrite = ',filewrite
      open(12,file=filewrite,form='unformatted',status='unknown')
      write(12)pix,pix1
      close(12)
*
      close(13)
*
*
      stop
      end
******
***********************************************************************
******
      subroutine gauss(isig,iRn,ia,ib,xxx,yyy,value,pix,ipix)
*
* modified: disk of constant brightness 
*
*
*
      integer isig,i1,i2,isig3sq,ix,iy,ipix,iRn,test,pixarea	
      real      pix(ipix,ipix)
      double precision normfac,factorex,sigsq2,value,dum,xxx,yyy
     &	      ,delx,dely,dist2,dist3,ia,ib
*
	
	ix    =  nint(xxx)
	iy    =  nint(yyy)
	delx  =  xxx - dble(ix)
	dely  =  yyy - dble(iy)
*
*
* 030711: isig3 = RADIUS  for source!!!
*
*
	isig3 = 3*isig 
c	isig3 =   isig 
	isig3sq = isig3**2
 	sigsq2  = 2.0*dble(isig)**2
	value = 0.0
	
	test=300
	pixarea=0	
c	 write(*,*)' in GAUSS'
cwrite(*,*)' isig,ix,iy,value = ',isig,ix,iy,value 
c	pause
*
* Gaussian source profiles with width  isig:
*
	      normfac = 0.0
	      do i1 = -isig3,isig3
		if(ix.eq.test)then
		write(*,*)
		endif 
	         do i2 = -isig3,isig3
*
c	            dist2 = i2*i2 + i1*i1
c		    dist3 = (i2-ib)*(i2-ib) + (i1-ia)*(i1-ia)		
*
	            dist2 = (i1-delx)**2  + (i2-dely)**2
		    dist3 = (i1-delx-ia)**2  + (i2-dely-ib)**2	
*
*
*
	            if(isig.ne.0)then
  	               if(dist2.le.dble(isig**2))then	
c 	                   factorex = exp(-dist2/sigsq2)
	                   factorex = 1.0
			  if(dist3.le.dble(iRn**2))then
			   factorex = 0.0
			  endif
			  if(factorex.eq.1.0)then
			   pixarea=pixarea+1
			  endif	 		
			   dum = pix(ix+i1,iy+i2)
	                   value = value + dum*factorex
	                   normfac = normfac + factorex	
			elseif(dist2.ge.dble(isig**2))then 	
			   factorex = 0.0		
                           dum = pix(ix+i1,iy+i2)
	                   value = value + dum*factorex
	                   normfac = normfac + factorex
	               endif
	            elseif(isig.eq.0)then
                      factorex = 1.0
                      dum = pix(ix+i1,iy+i2)
	              value = value + dum*factorex
	              normfac = normfac + factorex
		    endif
		if(ix.eq.test)then
		write(*,"(i1)",advance="no")  int(factorex)
		endif
	         enddo
	      enddo
	if(ix.eq.test)then
		write(*,*)
		write(*,*)
		write(*,*) 'Area in pixels =', pixarea
		endif 
	
*
 	      value = value/normfac 
		
	return
	end
