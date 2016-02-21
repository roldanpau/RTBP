# fixed points
set label "$p_2$" at first 5.626637805710596e+00, first 1.205112785642116e-01
set label "$p_3$" at first 6.015698455630493e+00, first -8.291411343795281e-08
set label "$z$" at first 5.712, first 0.062

#set out "invmfld_it23.tex"
#set term epslatex color
#set title "Invariant manifolds of $p_2$ and $p_3$"
set xlabel "$x$"
set ylabel "$p_x$"

#"stmfld_it2.res" w l lt 1 lc 2 t "stable", \
#"unstmfld_it3_neg.res" w l lt 1 lc 1 not, \

# plot also the points p2, p3, and z
plot \
"unstmfld_it2.res" w l lt 1 lc 1 t "unstable" , \
"stmfld_it3.res" w l lt 1 lc 2 t "stable" , \
"-" w p pt 3 pointsize 1 lc 3 not
 5.626637805710596e+00  1.205112785642116e-01
 6.015698455630493e+00 -8.291411343795281e-08
 5.712 0.062
E

#set term pop
#set out
#pause -1
