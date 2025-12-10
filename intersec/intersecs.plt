set term tikz
set out "intersecs.tex"

set key top left

set xlabel "$\\tJ$"
set ylabel "$x$"

plot \
	"intersecs_unst_br1.res" u 1:7 w l t "$z_1$", \
    "intersecs_st_br1.res" u 1:7 w l t "$z_2$", \
	"intersecs_unst_br2.res_nokinks" u 1:7 w l t "$z_3$", \
    "intersecs_st_br2.res_nokinks" u 1:7 w l t "$z_4$"

set term pop
set out
