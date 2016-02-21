unset label 

# fixed points
set label "p1" at first -3.041422643659351e-01, 7.992422910732644e-01
set label "p2" at first -3.041422643659351e-01, -7.992422910732644e-01
#set label "z" at first -0.087032, first 0.0

set out "invmfld_ell05_SEC2.eps"
set term post eps
#set title "Invariant manifolds of fixed point"
set xlabel "x"
set ylabel "px"

f(x)=0.0

plot [:] [-0.85:0.85]\
"unstmfld_ell05_SEC2.res" w l lt 1 lc 1 not, \
"unstmfld_ell05_p2_SEC2.res" w l lt 1 lc 1 not, \
"unstmfld_ell05_SEC2.res" u ($1):(-$2) w l lt 1 lc 2 not, \
"unstmfld_ell05_p2_SEC2.res" u ($1):(-$2) w l lt 1 lc 2 not, \
f(x) w l lt 3 lc 3 not

# "unstmfld_ell085.res" using ($1):(-$2) w l lt 1 lc 2 t "stable mfld", \

#"unstmfld_neg.res" w l lt 1 lc 1 not, \
#"stmfld.res" w l lt 1 lc 3 not, \
#"stmfld_neg.res" w l lt 1 lc 3 not, \
#"po.res" w p pt 3 pointsize 1 lc 2 not

set term pop
set out
