#set out "omega_j.eps"
#set term postscript enhanced monochrome eps

set xlabel "H"
set ylabel
#set key left top

plot \
"omega_neg.res" u 1:(-2.0*$2) t "omega1p" w l, \
"omega_neg.res" u 1:(-2.0*$3) t "omega1m" w l, \
"omega_neg.res" u 1:(-2.0*$4) t "omega2p" w l, \
"omega_neg.res" u 1:(-2.0*$5) t "omega2m" w l

set term pop
set out 
