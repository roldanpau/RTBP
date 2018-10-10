unset label 

# fixed points
set label "p" at first 0.2, 6.903739666022393e-01
set label "p" at first 2*pi-0.2, 6.903739666022393e-01
set label "z1" at first pi, 6.969766582609713e-01-0.0002
set label "z2" at first pi, 6.916448140296166e-01+0.0002

#set out "invmfld_H\-1.6614_SECg.eps"
#set term post eps

#set out "invmfld_ell05_SEC2.eps"
#set term post eps
#set title "Invariant manifolds of fixed point"
set xlabel "l"
set ylabel "L"

f(x)=pi

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

set arrow from pi, graph 0 to pi, graph 1 nohead ls 3

plot [0:2*pi] \
"unstmfld_H-1.6614_SECg_br1.res" u ($1<0 ? $1+2*pi : $1):2 w lp ls 1 not, \
"unstmfld_H-1.6614_SECg_br2.res" u ($1<0 ? $1+2*pi : $1):2 w l ls 1 not, \
"unstmfld_H-1.6614_SECg_br1.res" u ($1<0 ? -$1 : 2*pi-$1):2 w l ls 2 not, \
"unstmfld_H-1.6614_SECg_br2.res" u ($1<0 ? -$1 : 2*pi-$1):2 w l ls 2 not

#"unstmfld_neg.res" w l lt 1 lc 1 not, \
#"stmfld.res" w l lt 1 lc 3 not, \
#"stmfld_neg.res" w l lt 1 lc 3 not, \
#"po.res" w p pt 3 pointsize 1 lc 2 not

set term pop
set out
