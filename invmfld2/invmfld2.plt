unset label

# fixed points
set label "p0" at first -3.221507692925773e+00, first -3.757823611945582e-01
set label "p1" at first -5.116695151997712e+00, first -1.891957686237784e-01
set label "p2" at first -5.919861224972697e+00, first -5.874500915770870e-02
set label "p3" at first -5.919861224972697e+00, first 5.874500915770870e-02
set label "p4" at first -5.116695151997712e+00, first 1.891957686237784e-01
set label "p5" at first -3.221507692925773e+00, first 3.757823611945582e-01

#set out "invmfld2.eps"
#set term post eps
#set title "Invariant manifolds of fixed point"
set xlabel "x"
set ylabel "px"

plot \
"unstmfld.res" w l lt 1 lc 1 not, \
"unstmfld_neg.res" w l lt 1 lc 1 not, \
"stmfld.res" w l lt 1 lc 3 not, \
"stmfld_neg.res" w l lt 1 lc 3 not, \
"po.res" w p pt 3 pointsize 1 lc 2 not
#set term pop
#set out
