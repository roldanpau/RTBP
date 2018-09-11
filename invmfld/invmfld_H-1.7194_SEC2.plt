unset label 

# fixed points
#set label "p0" at first -6.389031137499175e-02, 2.621258562749938e+00
#set label "p1" at first -6.389031137499175e-02, -2.621258562749938e+00
#set label "z" at first -0.087032, first 0.0

set out "invmfld_H\-1.7194_SEC2.eps"
set term post enhanced color eps

#set title "Invariant manifolds of fixed point"
set xlabel "x"
set ylabel "px"

f(x)=0.0

#"< echo -7.115565698852543e-02 2.457642042989651e+00" w p pt 2 lc 7 not, \
#"< echo -7.115565698852543e-02 -2.457642042989651e+00" w p pt 2 lc 6 not, \
#"unstmfld_H-1.3594_p1.res" w l lc 7 not, \
#"unstmfld_H-1.3594_p1.res" u ($1):(-$2) w l lc 6 not, \
#"unstmfld_H-1.3594_p1_branch2.res" w d  lc 7 not, \
#"unstmfld_H-1.3594_p1_branch2.res" u ($1):(-$2) w d  lc 6 not, \
#"stmfld_H-1.3594_p1.res" w l  lc 6 not, \
#"stmfld_H-1.3594_p1.res" u ($1):(-$2) w l  lc 7 not, \
#"stmfld_H-1.3594_p1_branch2.res" w d  lc 6 not, \
#"stmfld_H-1.3594_p1_branch2.res" u ($1):(-$2) w d  lc 7 not, \

#Note: 
#unstmfld_H-1.7194_p1.res" u ($1):(-$2) = stmfld of p_2
#stmfld_H-1.7194_p1.res" u ($1):(-$2) = unstmfld of p_2

plot \
"< echo -4.470843776035918e-01 2.821703012713745e-01" w p pt 2 lc 7 not, \
"< echo -4.470843776035918e-01 -2.821703012713745e-01" w p pt 2 lc 6 not, \
"unstmfld_H-1.7194_p1.res" w l  lc 7 not, \
"unstmfld_H-1.7194_p1.res" u ($1):(-$2) w l  lc 6 not, \
"unstmfld_H-1.7194_p1_branch2.res" w d  lc 7 not, \
"unstmfld_H-1.7194_p1_branch2.res" u ($1):(-$2) w d  lc 6 not, \
"stmfld_H-1.7194_p1.res" w l  lc 6 not, \
"stmfld_H-1.7194_p1.res" u ($1):(-$2) w l  lc 7 not, \
"stmfld_H-1.7194_p1_branch2.res" w d  lc 6 not, \
"stmfld_H-1.7194_p1_branch2.res" u ($1):(-$2) w d  lc 7 not, \
f(x) w l lt 3 lc 2 not

#"unstmfld_neg.res" w l lt 1 lc 7 not, \
#"stmfld.res" w l lt 1 lc 6 not, \
#"stmfld_neg.res" w l lt 1 lc 6 not, \
#"po.res" w p pt 3 pointsize 1 lc 2 not

set term pop
set out
