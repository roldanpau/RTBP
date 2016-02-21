# fixed points
#set label "$p_0$" at first 4.204300492821935e+00, first 3.668601668676399e-01
#set label "$p_1$" at first 2.078884814357651e+00, first  3.668601668676399e-01

set xlabel "g"
set ylabel "G"

#plot "unstmfld_ell085.res" u 3:4 w l lc 1 t "unstable manifold", \
#"stmfld_ell085.res" u 3:4 w l lc 2 t "stable manifold"

plot \
"unstmfld_H-1.7194_p1.res" u 3:4 w l lc 1 t "unstable manifold", \
"stmfld_H-1.7194_p1.res" u 3:4 w l lc 3 t "stable manifold", \
"stmfld_H-1.7194_p2.res" u 3:4 w l lc 3 not, \
"unstmfld_H-1.7194_p2.res" u 3:4 w l lc 1 not

#"-" w p pt 3 pointsize 1 lc 3 not
#4.204300492821935e+00 3.668601668676399e-01
#2.078884814357651e+00 3.668601668676399e-01
#E
