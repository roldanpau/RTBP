set term tikz
set out "E2_bound.tex"

set xlabel "$\\tJ$"
set ylabel "$\\overline{E_2}$"

set key left top offset 5,0

plot \
"ebound_unst_br1.res" u 1:3 t "$i=1$" w l, \
"ebound_st_br1.res" u 1:3 t "$i=2$" w l, \
"ebound_unst_br2.res" u 1:3 t "$i=3$" w l, \
"ebound_st_br2.res" u 1:3 t "$i=4$" w l

set term pop
set out 
