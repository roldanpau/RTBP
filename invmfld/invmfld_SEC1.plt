#unset label 

# fixed points
#set label "p1" at first 2.190447719349280e-01, first 0.000000000000000e+00
#set label "p2" at first 7.216896649191896e-01, first 0
#set label "z" at first -0.087032, first 0.0

#set out "invmfld_ell05_SEC1.eps"
#set term post eps
#set title "Invariant manifolds of fixed point"
set xlabel "x"
set ylabel "px"

f(x)=0.0

plot \
"unstmfld_H-1.7194_SEC1.res" w lp lt 1 lc 1 t "unstable mfld", \
f(x) w l lt 3 lc 3 t "symmetry axis"

# "unstmfld_ell085.res" using ($1):(-$2) w l lt 1 lc 2 t "stable mfld", \

#"unstmfld_neg.res" w l lt 1 lc 1 not, \
#"stmfld.res" w l lt 1 lc 3 not, \
#"stmfld_neg.res" w l lt 1 lc 3 not, \
#"po.res" w p pt 3 pointsize 1 lc 2 not

#set term pop
#set out
