************************************************************************
************************************************************************
      PROGRAM MICROLENS
************************************************************************
************************************************************************
*
*  description of some integer variables:
*
*           lensmax     - maximal possible number of lenses
*           cellmax     - maximal possible number of cells
*           levelmax    - maximal possible number of levels
*           nlens       - number of lenses  (input)
*           ncell       - total (actual) number of cells (~ 1.7*nlens)
*           icell(i)    - number of lenses in cell with index i
*           indcell(i   - "actual" number of cell with index i
*           cellused    - number of cells used with cell-shooting
*           levelcel(i) - number of cells at level i with > 1 lens
*           levellen(i) - number of cells at level i with 1 lens
*           stack(i)    - index of cell at stack position i
*
*  description of some real variables:
*
*           xxx, yyy        - coordinates of ray in lens plane
*           xsidum, etadum  - coordinates of light ray in source plane
*           arand           - random seed number (input)
*           lensdata(1:3,i) -  values for x, y and mass of lens  #i     
*           cr1, cr2  - auxillary seed numbers for subroutine rands     
*           cell(1,i) - x-coordinate of c.m. of cell #i
*           cell(2,i) - y-coordinate of c.m. of cell #i
*           cell(3,i) - total mass of cell #i
*           cell(4,i) - sidelength of  cell #i
*           cell(5,i), cell(6,i)   - quadrupole moment of cell #i
*           cell(7,i), cell(8,i)   - octopole   moment of cell #i
*           cell(9,i), cell(10,i)  - 16-pole   moment of cell #i
*           cell(11,i), cell(12,i) - 32-pole   moment of cell #i
*           cell(13,i), cell(14,i) - 64-pole   moment of cell #i
*           distcomp2(i) - (sidelength of cells in level #i )**2
*           eps  - value of accuracy parameter
*           eps2 - eps**2
*
* use of subroutines:
*           input, setstar, numstar, setcell,
*           shootcl, shootfix, detko, shotayco
*	    output
*
* use of functions:
*           none
*					October 31, 1991; J. Wambsganss
*
*********************************************************************** 
      implicit none
      external input,setstar,numstar,setcell,shootcl,shootfix,detko
     &,shotayco,output
      external etime
*
      integer lensmax,cellmax,levelmax,ipix,ipix1,factor1,factor2
c     parameter (lensmax= 100000,cellmax=lensmax,levelmax=30
      parameter (lensmax= 10000,cellmax=lensmax,levelmax=15
     &,ipix=1000,ipix1=500,factor1=11,factor2=10)
*
      double precision xxx,yyy,xxxx,yyyy,del,xur,yur
     &,stamin,stadif,stafac,arand,arand0,cr1,cr2
     &,cell(14,cellmax),lensdata(3,lensmax),distcom2(levelmax),eps
     &,eps2,xsidum(factor2,factor2),etadum(factor2,factor2),raydist0
     &,facipix,pixdifh,xxxxx(factor2,factor2),yyyyy(factor2,factor2)
     &,sigma,gamma,agam,bgam,boa,bmsoams,rstars,sigmas,sigmac,raydist1
     &,pixmax0,pixdif,fracpixd,pixminx,pixmaxx,pixminy,raydist2
     &,pixmaxy,raydifx,rayminx,raymaxx,rayminy,raymaxy,raydist0sq
     &,avcell,avlens,alf(8),phi(8),alfx(factor1,factor1)
     &,alfy(factor1,factor1),xarr(factor1,factor1),yarr(factor1,factor1)
     &,alphaxar(factor1,factor1),alphayar(factor1,factor1),alphax,alphay
     &,masstot,massaver,minmass,maxmass,power,powerp1,powerp2,powerp3
     &,ampth,pixlen(2),aaa1,aaa2
*
      integer nlens,ncell,nray,nrayx,nrayy,ix,iy,indlens(lensmax)
     &,levelcel(levelmax),levellen(levelmax),stack(cellmax)
     &,highle,debug,clno(lensmax),totused,lno(lensmax),jxmin,jxmax,jymax
     &,jymin,jobnu,        ix1,iy1,icell(cellmax),cellused
     &,jxdum(factor2*factor2),jydum(factor2*factor2),indcell(cellmax)
     &,ixdum(factor2*factor2),iydum(factor2*factor2),lensused
      integer*8 rayshot 
*
      integer*4 pix(ipix,ipix),pix1(ipix1,ipix1),month,day,year
     &       ,month1,day1,year1
*
      integer*2 pixa(ipix,ipix),pixb(ipix1,ipix1)
*
      real   cpustart, cpufinish
      real time1,time2,t1,t2,times(20)
*
      character*3 string1
      character*8 string
      character*2 string_vec(3) 
      integer*4 ivec1(3),ivec2(3)
      integer values(8)
*
* COMMON-blocks: 
*    "test" in each subroutine and function; 
*    "m"    in main, massf, massav;
*    "leng" in main,input,numstar,numcell; 
*    "pix"  in main, raycount;
*    "time" in main,  output, setstar, setcell, shotayco
*
      common/test/debug
      common/m/minmass,maxmass,power,powerp1,powerp2,powerp3,aaa1,aaa2
      common/leng/stamin,stadif
      common/pix/facipix,pixdifh
      common/jkwtime/times
*
* determine day, month, year;  system time (CONVEX routines!):  
*
*
	call date_and_time(VALUES=values)
	day   = values(3)
	month = values(2)
	year  = values(1)
*
	write(string_vec(1),'(i2.2)')values(5)
	write(string_vec(2),'(i2.2)')values(6)
	write(string_vec(3),'(i2.2)')values(7)
	string = string_vec(1)//':'//string_vec(2)//':'//string_vec(3)
*
      call input(debug,nray,arand,sigmas,sigmac,gamma,eps,minmass
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
*
* determine coordinates (and masses) of all lenses:                 
*
      call cpu_time ( t1 )
      cpustart  = t1
      if(nlens.gt.0)then
         call setstar(arand,cr1,cr2,rstars,nlens,lensmax,lensdata
     &			,masstot,pixdifh)
      endif
      call cpu_time ( t2 )
      times(1)  = t2-t1
      t1        = t2
*
* assign numbers to each lens:
*
      if(nlens.gt.0)then
         call numstar(nlens,lensmax,lensdata,levelmax,indlens)          
      endif
      call cpu_time ( t2 )
      times(2)  = t2-t1
      t1        = t2
*
* build these lenses into cell structure:                           
*
      if(nlens.gt.0)then
         call setcell(cellmax,lensmax,levelmax,cell,icell,indcell
     &	,ncell,distcom2,lensdata,indlens,nlens,highle,levelcel,levellen)
      endif
      call cpu_time ( t2 )
      times(3) = t2-t1
      t1       = t2 
*
* start shooting:
*
      write(*,*)' shooting starts ...'
      t1    = 0.0d0
      t2    = 0.0d0
      time1 = 0.0d0
      time2 = 0.0d0
      call cpu_time ( time1 )
*
*********************** RAY SHOOTING STARTS ****************************
*
      do 100  ix = jxmin,jxmax
         xxx    =    raydist0 * float(ix)
         do 200  iy  = jymin, jymax
            yyy    =   raydist0 * float(iy)
*
	    if(iy.eq.jymin)then
      call cpu_time ( cpustart )
        write(*,1095)ix,iy,xxx,yyy,cpustart
 1095 format(' ix,iy,xxx,yyy,CPUsec:',2i6,2f15.5,f15.1)
 1096 format('   ix,iy,xxx,yyy,CPUsec:',2i6,2f15.5,f15.1)
	    endif
*
* determination of cell/lens structure for the bunch of rays centered
* on (xxx/yyy), no actual shooting is done:
*
      call cpu_time ( t1 )
            call shootcl(cellmax,lensmax,stack,xxx,yyy,cell,eps2
     &         ,cellused,icell,indcell,ncell,nlens,indlens,lensdata
     &         ,lensused,clno,lno,raydist0sq,totused)
      call cpu_time ( t2 )
      	    times(4)  = times(4) + (t2-t1)
*
* shoot test rays in level 2 (only fixed cells, NO lenses,
* without shear (!), to determine Taylor coefficients)
*
            do ix1= 1,factor1
               do iy1= 1,factor1                        
                  xxxx = xxx - 0.5*raydist0 + (ix1-1)*raydist1
                  yyyy = yyy - 0.5*raydist0 + (iy1-1)*raydist1
*
      		call cpu_time ( t1 )
                  call shootfix(cellmax,lensmax,xxxx,yyyy,alphax
     &	             ,alphay,cell,lensdata, lno,lensused,totused)
      		call cpu_time ( t2 )
      		  times(5)  = times(5) + (t2-t1)
*
                  xarr(ix1,iy1) = xxxx 
                  yarr(ix1,iy1) = yyyy
                  alfx(ix1,iy1) = alphax
                  alfy(ix1,iy1) = alphay
               enddo
	    enddo
*
* deflection angle due to cells at four corners:
*
            do ix1= 1,factor1-1                           
               do iy1= 1,factor1-1                      
                  alf(1) = alfx(ix1+1,iy1+1)               
                  alf(2) = alfy(ix1+1,iy1+1)               
                  alf(3) = alfx(ix1  ,iy1+1)               
                  alf(4) = alfy(ix1  ,iy1+1)               
                  alf(5) = alfx(ix1  ,iy1  )               
                  alf(6) = alfy(ix1  ,iy1  )               
                  alf(7) = alfx(ix1+1,iy1  )               
                  alf(8) = alfy(ix1+1,iy1  )               
                  xxxx   = xarr(ix1,iy1)                   
                  yyyy   = yarr(ix1,iy1)                   
                  call detko(alf,phi,del)                  
      		call cpu_time ( t1 )
*
* shooting in level 3: lenses directly included, cell contribution 
*  interpolated:
*
                  call shotayco(phi,del,lensmax,lensdata,lno,xxxx,yyyy
     &		   ,agam,bgam,sigmac,sigmas,lensused,pix,ipix,pixminx
     &             ,pixmaxx,pixminy,pixmaxy,pix1,ipix1,factor2,rayshot
     &		   ,raydist2,xsidum,etadum,xxxxx,yyyyy,alphaxar,alphayar
     &             ,ixdum,iydum,jxdum,jydum,nlens)
      		call cpu_time ( t2 )
                 times(6) = times(6) + (t2-t1)
               enddo
            enddo
*
            avcell = avcell + cellused
            avlens = avlens + lensused
*
  200    continue
  100 continue
*
********************* RAY SHOOTING ENDS ********************************
*
      write(*,*)' ..............shooting ends'
*
      		call cpu_time ( time2 )
      times(7)   = times(7) + time2 - time1
*
      times(8) = time2
*
      call output(ampth,rayshot,avlens,avcell,levelmax,levelcel
     &,levellen,highle,nlens,rstars,ncell,xur,yur,lensdata,lensmax
     &,indlens,cell,cellmax,indcell,pixlen,raydist0,raydist1,raydist2
     &,nray,nrayx,nrayy,rayminx,raymaxx,rayminy,raymaxy,jxmin,jxmax
     &,jymin,jymax,raydifx,boa,bmsoams,factor1,factor2,jobnu,arand0
     &,debug,sigmas,sigmac,gamma,eps,minmass, maxmass,power,massaver
     &,masstot,pixmax0,pixminx,pixminy,pixdif,fracpixd,pix,pix1,ipix
     &,ipix1,pixa,pixb,month,day,year,string,string1)
*
      stop
      end
