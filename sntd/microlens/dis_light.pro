; 
;       idl program for visualization of microlens results
;		Joachim Wambsganss, January 15, 2008
;
device, decompose = 0
loadct, 4
;
;
a = intarr(1000,1000,/nozero)
b = intarr(500,500,/nozero)
;
;
;openr, 1, 'IRIS401'
openr, 1, 'IRIS401-track'
readu, 1, a,b
close, 1
;
;
;
window, 1, xpos=  70, ypos= 50, title='Microlens: Magnification Pattern', xsize=1000, ysize=1000
tv, BYTSCL(a,min=900,max=1750)
;tv, BYTSCL(a)
;
window, 3, xpos=  1450, ypos= 750, title='Microlens: Check', xsize=500, ysize=500
tv, BYTSCL(b)
tv, BYTSCL(b,min=700,max=1750)
;



nlines=FILE_LINES('out_line')
openr,lun,"out_line",/GET_LUN
dum1     = 0.0
dum2     = 0.0
t        = 0.0
xval     = 0.0
yval     = 0.0
pixvalue = 0.0
maglin   = 0.0
xpix     = 0.0
ypix     = 0.0
time=fltarr(nlines)
mag=fltarr(nlines)
	FOR j=0,nlines-1 DO BEGIN
	readf,lun, t,xval,yval,pixvalue,maglin,xpix,ypix
	time(j)=t
	mag(j)=maglin
	ENDFOR
close,lun
;
dum1=time(0)
dum2=time(nlines-1)
window, 4, xpos=   70, ypos= 950, title='Lightcurve    ', xsize=1000, ysize=300
plot,time,mag,XRANGE=[dum1,dum2],YRANGE=[0,6],TITLE='light curve',xtitle='pixels ',ytitle='magnification [linear]',/xstyle


;
end
