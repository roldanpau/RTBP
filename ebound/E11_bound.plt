set term tikz
set out "E11_bound.tex"

set xlabel "$\\tJ$"
set ylabel "$\\overline{\\EE_1^1}$"

set key left top offset 5,0

plot \
"ebound_unst_br1.res" u 1:2 t "$i=1$" w l, \
"ebound_st_br1.res" u 1:2 t "$i=2$" w l, \
"ebound_unst_br2.res" u 1:2 t "$i=3$" w l, \
"ebound_st_br2.res" u 1:2 t "$i=4$" w l

set term pop
set out 
