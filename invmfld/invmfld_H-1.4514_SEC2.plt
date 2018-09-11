unset label 

# fixed points
#set label "p0" at first -6.389031137499175e-02, 2.621258562749938e+00
#set label "p1" at first -6.389031137499175e-02, -2.621258562749938e+00
#set label "z" at first -0.087032, first 0.0

#set out "invmfld_H\-1.4514_SEC2.eps"
#set term post eps

#set out "invmfld_ell05_SEC2.eps"
#set term post eps
#set title "Invariant manifolds of fixed point"
set xlabel "x"
set ylabel "px"

f(x)=0.0

#"< echo -7.115565698852543e-02 2.457642042989651e+00" w p pt 2 lc 1 not, \
#"< echo -7.115565698852543e-02 -2.457642042989651e+00" w p pt 2 lc 3 not, \
#"unstmfld_H-1.3594_p1.res" w l lc 1 not, \
#"unstmfld_H-1.3594_p1.res" u ($1):(-$2) w l lc 3 not, \
#"unstmfld_H-1.3594_p1_branch2.res" w d  lc 1 not, \
#"unstmfld_H-1.3594_p1_branch2.res" u ($1):(-$2) w d  lc 3 not, \
#"stmfld_H-1.3594_p1.res" w l  lc 3 not, \
#"stmfld_H-1.3594_p1.res" u ($1):(-$2) w l  lc 1 not, \
#"stmfld_H-1.3594_p1_branch2.res" w d  lc 3 not, \
#"stmfld_H-1.3594_p1_branch2.res" u ($1):(-$2) w d  lc 1 not, \

#Note: 
#unstmfld_H-1.4514_p1.res" u ($1):(-$2) = stmfld of p_2
#stmfld_H-1.4514_p1.res" u ($1):(-$2) = unstmfld of p_2

plot \
"unstmfld_H-1.4514_p1.res" w lp  lc 1 not, \
"unstmfld_H-1.4514_p1_branch2.res" w d  lc 1 not, \
f(x) w l lt 3 lc 2 not

#"unstmfld_H-1.4514_p1.res" u ($1):(-$2) w l  lc 3 not, \
#"unstmfld_H-1.4514_p1_branch2.res" u ($1):(-$2) w d  lc 3 not, \
#"stmfld_H-1.4514_p1.res" w l  lc 3 not, \
#"stmfld_H-1.4514_p1.res" u ($1):(-$2) w l  lc 1 not, \
#"stmfld_H-1.4514_p1_branch2.res" w d  lc 3 not, \
#"stmfld_H-1.4514_p1_branch2.res" u ($1):(-$2) w d  lc 1 not, \

#"unstmfld_neg.res" w l lt 1 lc 1 not, \
#"stmfld.res" w l lt 1 lc 3 not, \
#"stmfld_neg.res" w l lt 1 lc 3 not, \
#"po.res" w p pt 3 pointsize 1 lc 2 not

set term pop
set out
