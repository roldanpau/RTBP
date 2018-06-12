set out "omega_j_br2.eps"
set term postscript enhanced monochrome eps

set xlabel "H"
set ylabel
#set key left top

plot "omega_neg.res" u 1:(-2.0*$3) t "omega1" w l, \
"omega_neg.res" u 1:(-2.0*$5) t "omega2" w l

set term pop
set out 
