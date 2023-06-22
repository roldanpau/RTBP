#set out "omega_j.eps"
#set term postscript enhanced color eps

set xlabel "H"
set ylabel "alpha"

#set key left top
set key spacing 3

set datafile missing "?"

plot \
"omega_neg.res" u 1:(-2.0*$2) t "alphau1" w l, \
"omega_neg.res" u 1:(-2.0*$3) t "alphau2" w l, \
"omega_neg.res" u 1:(2.0*$4) t "alphas1" w l, \
"omega_neg.res" u 1:(2.0*$5) t "alphas2" w l

set term pop
set out 
