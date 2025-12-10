set term tikz
set out "split.tex"

set xlabel "$\\tJ$"
set ylabel "$\\theta_i$"

set key top right

#f(x)=0.0

set grid

plot \
"splittings_unst_br1.res" w l t "$\\splitting_1$", \
"splittings_st_br1.res" w l t "$\\splitting_2$", \
"splittings_unst_br2.res" w l t "$\\splitting_3$", \
"splittings_st_br2.res" w l t "$\\splitting_4$"

#f(x) w l lt 3 lc 3 not

set term pop
set out
