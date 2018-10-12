#set out "omega_neg.eps"
#set term eps color
set xlabel "H"
set ylabel "omega_neg"
set key left top
plot "omega_neg.res" u 1:2 t "omega_neg_SECg_br1" w l, \
"omega_neg.res" u 1:3 t "omega_neg_SECg_br2" w l, \
"omega_neg.res" u 1:4 t "omega_neg_SECg2_br1" w l, \
"omega_neg.res" u 1:5 t "omega_neg_SECg2_br2" w l
set term pop
set out 
