unset label

# fixed points
set label "p2" at first -5.919861224972697e+00, first -5.874500915770870e-02
set label "p3" at first -5.919861224972697e+00, first 5.874500915770870e-02
set label "Ws2p3" at first -5.803, first 0.0494
set label "Wu1p3" at first -6.148, first 0.0494
set label "Ws1p2" at first -6.148, first -0.0498
set label "Wu2p2" at first -5.807, first -0.0498
#set label "$z$" at first -6.211, first 0.0
#set label "$\\hat z$" at first -5.813, first 0.0

#set out "invmfld2_it23.eps"
#set term post eps
#set title "Invariant manifolds of $p_2$ and $p_3$"
set xlabel "x"
set ylabel "px"

# plot also the points p2, p3, and intersection points
plot [-6.5:-5.5] [-0.07:0.07] \
"unstmfld.res" w l lt 1 lc 1 not, \
"unstmfld_neg.res" w l lt 1 lc 1 not, \
"stmfld.res" w l lt 1 lc 3 not, \
"stmfld_neg.res" w l lt 1 lc 3 not, \
"-" w p pt 3 pointsize 1 lc 2 not
-5.919861224972697e+00 -5.874500915770870e-02
-5.919861224972697e+00 5.874500915770870e-02
-6.211 0.0
-6.161 0.0
-6.126 0.0
-5.902 0.0
-5.861 0.0
-5.813 0.0
E

#set term pop
#set out

