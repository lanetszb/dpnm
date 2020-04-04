# set term png  
# set output "fig_press_flowrates.png"
set term tikz latex size 9,6
set output 'fig3.txt'

set xlabel "Time, $t$ ($s$)"
set ylabel "Average pore pressure, $p$ ($10^{-3} Pa$)" offset 1.7,0,0
set y2label "Flow rate, $Q$ ($10^{3} \\; m^3/sec$)" offset +1.3,0,0
set autoscale
# set yrange [299000:299300]
# set ytics 0.5 offset .9,0,0
#set y2range [1.018e-5:1.028e-5]
# set y2tics 5 offset -.9,0,0
# set xrange [0:15]
set xtics 3 offset 0,.1,0 
unset grid
set xtics nomirror
set ytics nomirror
set y2tics

plot 'fig_press_flowrates.txt' using 1:2 with l lc 1 title 'Average pore pressure', \
'' using 1:3 with l lc 3 axes x1y2  title 'Inlet flow rate', \
'' using 1:4 with l lc 4 axes x1y2  title 'Outlet flow rate'