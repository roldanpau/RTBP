#set term tikz
#set out "Lbound.tex"

set xlabel "$\\tJ$"
set ylabel "$\\overline{L_i}$"

set key left top spacing 2

plot [:-1.485] \
"Lbound_unst_br1.res" u 1:2 t "$\\overline{L_1}$" w l, \
"Lbound_st_br1.res" u 1:2 t "$\\overline{L_2}$" w l, \
"Lbound_unst_br2.res" u 1:2 t "$\\overline{L_3}$" w l, \
"Lbound_st_br2.res" u 1:2 t "$\\overline{L_4}$" w l

#plot [:-1.485] \
#"Lbound_unst_br1.res" u 1:3 t "$\\overline{a_1}$" w l, \
#"Lbound_st_br1.res" u 1:3 t "$\\overline{a_2}$" w l, \
#"Lbound_unst_br2.res" u 1:3 t "$\\overline{a_3}$" w l, \
#"Lbound_st_br2.res" u 1:3 t "$\\overline{a_4}$" w l

set term pop
set out 
