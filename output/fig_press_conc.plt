# set term png  
# set output "fig_press_conc.png"

set term tikz latex size 9,6
set output 'fig2.txt'
set xlabel "Time, $t$ ($s$)"
set ylabel "Average pore pressure, $p$ ($10^{-3} Pa$)" offset 1.7,0,0
set y2label "Average matrix concentration, $C_{m}$ ($10^{3} \\; kg/m^3$)" offset +1.3,0,0
set autoscale
# set yrange [299000:299300]
# set autoscale y
# set ytics 0.5 offset .9,0,0
set y2range [1.45:1.5]
# set y2tics 5 offset -.9,0,0
# set xrange [0:15]
set xtics 3 offset 0,.1,0 
unset grid
set xtics nomirror
set ytics nomirror
set y2tics

plot 'fig_press_conc.txt' using 1:2 with l lc 1 title 'Average pore pressure', \
'' using 1:3 with l lc 3 axes x1y2  title 'Average matrix concentration'