#set out "omega_j_zoom2.eps"
#set term postscript enhanced color eps

set xlabel "H"
set ylabel
#set key left top
set key spacing 3

# To zoom:

plot [-1.6125:-1.59]\
"omega_neg.res" u 1:(-2.0*$2) t "omega1p" w lp, \
"omega_neg.res" u 1:(-2.0*$3) t "omega1m" w lp, \
"omega_neg.res" u 1:(-2.0*$4) t "omega2p" w lp, \
"omega_neg.res" u 1:(-2.0*$5) t "omega2m" w lp

set term pop
set out 
