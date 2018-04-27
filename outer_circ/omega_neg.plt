set out "omega_neg.eps"
set term eps color
set xlabel "H"
set ylabel "omega_neg"
plot "omega_neg_f_SECg.res" u 1:2 t "omega_neg^1" w lp, \
"omega_neg_f_SECg.res" u 1:3 t "omega_neg^2" w lp
set term pop
set out 
