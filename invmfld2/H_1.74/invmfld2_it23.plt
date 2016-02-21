# fixed points
set label "p2" at first -5.560849822621776e+00, first -5.681451293281249e-02
set label "p3" at first -5.560849822621776e+00, first 5.681451293281249e-02
set label "z1" at first -5.73, first 0.0
set label "z2" at first -5.58, first 0.0
set label "Ws2p3" at first -5.521,0.0405
set label "Wu1p3" at first -5.664,0.0405
set label "Ws1p2" at first -5.664,-0.0405
set label "Wu2p2" at first -5.521,-0.0405

#set out "invmfld2_it23_H174.eps"
#set term post eps
#set title "Invariant manifolds of $p_2$ and $p_3$"
set xlabel "x"
set ylabel "px"

# plot also the points p2, p3
plot [-5.8:-5.45] [-0.07:0.07] \
"unstmfld.res" w l lt 1 lc 1 not, \
"unstmfld_neg.res" w l lt 1 lc 1 not, \
"stmfld.res" w l lt 1 lc 3 not, \
"stmfld_neg.res" w l lt 1 lc 3 not , \
"-" w p pt 3 pointsize 1 lc 2 not
-5.560849822621776e+00 -5.681451293281249e-02
-5.560849822621776e+00 5.681451293281249e-02
E

#set term pop
#set out

