#figure 'I4D6z122_Re.ps' 
set style data linespoints
set xlabel "x" 
#set ylabel " " 
set logscale x
#set title 'Re I_4^(n+2)(z1*z2^2)' 
#set yr [-.003:-.0018] 
plot \
"demo_detG.dat" using 1:2 title "Re I_4^(n+2)(z1*z2^2)" 
# "demo_detG.dat" using 1:3 title "Im I_4^(n+2)(z1*z2^2)" 
set terminal postscript color
set output 'I4D6z122_Re.ps'
replot
set output
set terminal x11
