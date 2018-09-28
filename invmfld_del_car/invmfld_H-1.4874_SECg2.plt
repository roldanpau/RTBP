unset label 

# fixed points
set label "p" at first -pi+0.2, 6.923540823087978e-01
set label "p" at first pi-0.2, 6.923540823087978e-01
set label "z1" at first 0, 6.977129268113177e-01-0.0002
set label "z2" at first 0, 6.928500635543233e-01+0.0002

#set out "invmfld_H\-1.4874_SECg2.eps"
#set term post eps
#set title "Invariant manifolds of fixed point"
set xlabel "l"
set ylabel "L"

f(x)=0.0

set style line 1 lt 10 lc rgb "red" 
set style line 2 lt 0 lc rgb "blue" 
set style line 3 lt 5 lc rgb "gray" 

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
#unstmfld_H-1.7194_p1.res" u ($1):(-$2) = stmfld of p_2
#stmfld_H-1.7194_p1.res" u ($1):(-$2) = unstmfld of p_2

set arrow from 0, graph 0 to 0, graph 1 nohead ls 3

#plot [-pi:pi] \

plot [:] \
"unstmfld_H-1.4874_SECg2_br1.res" u 1:2 w l ls 1 not, \
"unstmfld_H-1.4874_SECg2_br2.res" u 1:2 w lp ls 1 not, \
"unstmfld_H-1.4874_SECg2_br1.res" u (-$1):2 w l ls 2 not, \
"unstmfld_H-1.4874_SECg2_br2.res" u (-$1):2 w lp ls 2 not

#plot [-0.1:2*pi+0.1] \
#"unstmfld_H-1.7194_SECg2_br1.res" u ($1<0 ? $1+2*pi : $1):2 w l ls 1 not, \
#"unstmfld_H-1.7194_SECg2_br2.res" u ($1<0 ? $1+2*pi : $1):2 w l ls 1 not, \
#"unstmfld_H-1.7194_SECg2_br1.res" u ($1<0 ? -$1 : 2*pi-$1):2 w l ls 2 not, \
#"unstmfld_H-1.7194_SECg2_br2.res" u ($1<0 ? -$1 : 2*pi-$1):2 w l ls 2 not

#"po.res" w p pt 3 pointsize 1 lc 2 not

set term pop
set out
