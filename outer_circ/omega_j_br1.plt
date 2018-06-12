set out "omega_j_br1.eps"
set term postscript enhanced monochrome eps

set xlabel "H"
set ylabel
#set key left top

plot "omega_neg.res" u 1:(-2.0*$2) t "omega1" w l, \
"omega_neg.res" u 1:(-2.0*$4) t "omega2" w l

set term pop
set out 
