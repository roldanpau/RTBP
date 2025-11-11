set out "E11_bound.eps"
set term postscript enhanced color eps

set xlabel "J"
set ylabel "E11"

#set key left top
#set key spacing 3

plot \
"ebound_unst_br1.res" u 1:2 t "u1" w l, \
"ebound_unst_br2.res" u 1:2 t "u2" w l, \
"ebound_st_br1.res" u 1:2 t "s1" w l, \
"ebound_st_br2.res" u 1:2 t "s2" w l

set term pop
set out 
