# set term png  
# set output "fig_langmuir.png"

set term tikz latex size 9,6
set output 'fig4.txt'
set xlabel "$P$, $Pa$"
set ylabel "$C$, $kg/m^3$" offset 0.01,0,0


a = 5.07024e-01
b = 1.43614e-06
c = -2.66609e-14
d = -2.05746e-20
e = 1.78042e-27

set xrange[0:300000]

f(x) = a*x**0 + b*x**1 + c*x**2 + d*x**3 + e*x**4

set ytics 0.1

unset grid
set xtics nomirror
set ytics nomirror


plot f(x) with lines lw 2 lc rgb "black" notitle