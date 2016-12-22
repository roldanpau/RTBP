unset label 

# fixed points
#set label "p0" at first -6.389031137499175e-02, 2.621258562749938e+00
#set label "p1" at first -6.389031137499175e-02, -2.621258562749938e+00
#set label "z" at first -0.087032, first 0.0

set out "invmfld_H\-1.3594_SEC2.eps"
set term post eps

#set out "invmfld_ell05_SEC2.eps"
#set term post eps
#set title "Invariant manifolds of fixed point"
set xlabel "g"
set ylabel "G"

f(x)=0.0

set style line 1 lt 10 lc rgb "red" 
set style line 2 lt 0 lc rgb "blue" 

plot [-pi-0.1:pi+0.1] \
"unstmfld_H-1.3594_p1_it1.res" w l ls 1 not, \
"unstmfld_H-1.3594_p1_it1.res" u (-$1):($2) w l ls 2 not, \
"unstmfld_H-1.3594_p1_it2.res" u ($1<pi ? $1 : 1/0):2 w l ls 1 not, \
"unstmfld_H-1.3594_p1_it2.res" u ($1<pi ? -$1 : 1/0):2 w l ls 2 not, \
"unstmfld_H-1.3594_p1_it3.res" u ($1-2*pi):($2) w l ls 1 not, \
"unstmfld_H-1.3594_p1_it3.res" u (2*pi-$1):($2) w l ls 2 not, \
"unstmfld_H-1.3594_p1_branch2_it1.res" u ($1-2*pi):($2) w l ls 1 not, \
"unstmfld_H-1.3594_p1_branch2_it1.res" u (2*pi-$1):($2) w l ls 2 not, \
"unstmfld_H-1.3594_p1_branch2_it2.res" w l ls 1 not, \
"unstmfld_H-1.3594_p1_branch2_it2.res" u (-$1):($2) w l ls 2 not, \
"unstmfld_H-1.3594_p1_branch2_it3.res" u ($1>pi ? $1-2*pi : 1/0):2 w l ls 1 not, \
"unstmfld_H-1.3594_p1_branch2_it3.res" u ($1>pi ? 2*pi-$1 : 1/0):2 w l ls 2 not

#plot [-pi-0.1:pi+0.1] \
#"< echo 1.045921189040063e+00 6.805459985644599e-01" w p pt 2 lc 0 not, \
#"< echo -1.045921189040063e+00 6.805459985644599e-01" w p pt 2 lc 0 not, \
#"< echo 3.141592655269312e+00 6.786291606785653e-01" w p pt 2 lc 0 not, \
#"< echo -3.141592655269312e+00 6.786291606785653e-01" w p pt 2 lc 0 not, \
#"unstmfld_H-1.3594_p1_it1.res" w l ls 1 not, \
#"unstmfld_H-1.3594_p1_it1.res" u (-$1):($2) w l ls 2 not, \
#"unstmfld_H-1.3594_p1_it2.res" w l ls 1 not, \
#"unstmfld_H-1.3594_p1_it2.res" u (-$1):($2) w l ls 2 not, \
#"unstmfld_H-1.3594_p1_it3.res" u ($1-2*pi):($2) w l ls 1 not, \
#"unstmfld_H-1.3594_p1_it3.res" u (2*pi-$1):($2) w l ls 2 not, \
#"unstmfld_H-1.3594_p1_branch2_it1.res" u ($1-2*pi):($2) w l ls 1 not, \
#"unstmfld_H-1.3594_p1_branch2_it1.res" u (2*pi-$1):($2) w l ls 2 not, \
#"unstmfld_H-1.3594_p1_branch2_it2.res" w l ls 1 not, \
#"unstmfld_H-1.3594_p1_branch2_it2.res" u (-$1):($2) w l ls 2 not, \
#"unstmfld_H-1.3594_p1_branch2_it3.res" u ($1>pi ? $1-2*pi : 1/0):2 w l ls 1 not, \
#"unstmfld_H-1.3594_p1_branch2_it3.res" u ($1>pi ? 2*pi-$1 : 1/0):2 w l ls 2 not

#"unstmfld_neg.res" w l lt 1 lc 1 not, \
#"stmfld.res" w l lt 1 lc 3 not, \
#"stmfld_neg.res" w l lt 1 lc 3 not, \
#"po.res" w p pt 3 pointsize 1 lc 2 not

set term pop
set out
