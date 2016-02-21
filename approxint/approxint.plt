# fixed points
set label "$p_2$" at first 5.626637805710596e+00, first 1.205112785642116e-01
set label "$z$" at first 5.777264327513651e+00, first 1.012387819028703e-05

#set title "Approximate intersection of W^u(p_2) with x axis"
set xlabel "$x$"
set ylabel "$p_x$"

set out "approxint.tex"
set term epslatex color

# Plot: 
# 1) unst manifold of p_2
# 2) x axis 
# 3) point p_2
# 4) approx intersection point z

plot [5.5:6] [-0.05:0.15] \
"unstmfld_it2.res" w l lc 1 t "unstable", \
0 w l lc 2 not, \
"-" w p pt 3 pointsize 1 lc 3 not
5.626637805710596e+00 1.205112785642116e-01
5.777264327513651e+00 1.012387819028703e-05
E

set term pop
set out
#pause -1
