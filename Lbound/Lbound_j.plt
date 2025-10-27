#set out "Lbound_j.eps"
#set term postscript enhanced color eps

set xlabel "H"
set ylabel "Lbound"

#set key left top
#set key spacing 3

plot \
"Lbound_unst_br1.res" u 1:2 t "u1" w l, \
"Lbound_unst_br2.res" u 1:2 t "u2" w l, \
"Lbound_st_br1.res" u 1:2 t "s1" w l, \
"Lbound_st_br2.res" u 1:2 t "s2" w l

set term pop
set out 
