#figure 'demo_a55_dets_sing.ps' 
set style data linespoints
set xlabel "pt5" 
#set ylabel " " 
set logscale x
#set title 'scattering singularity' 
#set yr [-.001:.001] 
plot \
"demo_a55_dets_sing.dat" using 1:3 title "Re[A55(1,1,1,1,1)] " ,\
"demo_a55_dets_sing.dat" using 1:4 title "Im[A55(1,1,1,1,1)] " 
set terminal postscript color
set output 'demo_a55_dets_sing.ps'
replot
set output
set terminal x11
