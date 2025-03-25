#figure 'I4D6z122_Im.ps' 
set style data linespoints
set xlabel "x" 
#set ylabel " " 
set logscale x
#set title 'Im I_4^(n+2)(z1*z2^2)' 
#set xr [xmin:xmax] 
#set yr [ymin:ymax] 
plot \
 "demo_detG.dat" using 1:3 title "Im I_4^(n+2)(z1*z2^2)" 
set terminal postscript color
set output 'I4D6z122_Im.ps'
replot
set output
set terminal x11
