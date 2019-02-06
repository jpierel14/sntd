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
openr, 1, 'IRIS401'
readu, 1, a,b
close, 1
;
;
;
window, 1, xpos=  70, ypos= 50, title='Microlens: Magnification Pattern', xsize=1000, ysize=1000
tv, BYTSCL(a,min=800,max=1750)
;tv, BYTSCL(a)
;
window, 3, xpos=  1450, ypos= 750, title='Microlens: Check', xsize=500, ysize=500
tv, BYTSCL(b)
;
;
end
