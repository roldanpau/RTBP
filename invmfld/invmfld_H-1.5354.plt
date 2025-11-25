unset label 

# fixed points
#set label "p0" at first 1.423330969133599e-01,0.000000000000000e+00 \
#point pt 3 ps 2 offset 1,1
set label "p1" at first -1.869996395538581e-01,1.266386459013741e+00 \
point pt 3 ps 2 offset 1,1
#set label "p2" at first 8.186568343573867e-01,1.354558826216490e-14 \
#point pt 3 ps 2 offset -4,1
set label "p3" at first -1.869996395538581e-01,-1.266386459013746e+00 \
point pt 3 ps 2 offset 1,-1
#set label "p1" at first -6.389031137499175e-02, -2.621258562749938e+00
#set label "z" at first -0.087032, first 0.0

set out "invmfld_H\-1.5354.eps"
set term post eps color

#set title "Invariant manifolds of fixed point"
set xlabel "x"
set ylabel "px"

f(x)=0.0

set style line 1 lt 10 lc rgb "red"
set style line 2 lt 0 lc rgb "blue" lw 2

#Note: 
#unstmfld_H-1.5354_p1.res" u ($1):(-$2) = stmfld of p_2
#stmfld_H-1.5354_p1.res" u ($1):(-$2) = unstmfld of p_2
# stmfld of p_1 branch1 = sym_rev(unstmfld of p_1 branch2)
# stmfld of p_1 branch2 = sym_rev(unstmfld of p_1 branch1)

#plot [-1:0]\

#"po_H-1.5354_p1.res" w p pt 2 ps 2 lc 0 not, \
#"po_H-1.5354_p2.res" w p pt 2 ps 2 lc 0 not, \

plot \
"unstmfld_H-1.5354_p1_branch1.res" u ($1):($3) w l ls 1 not, \
"unstmfld_H-1.5354_p1_branch2.res" u ($1):($3) w l ls 1 not, \
"stmfld_H-1.5354_p1_branch1.res" u ($1):($3) w l ls 2 not, \
"stmfld_H-1.5354_p1_branch2.res" u ($1):($3) w l ls 2 not, \
"unstmfld_H-1.5354_p1_branch1.res" u ($1):(-$3) w l ls 2 not, \
"unstmfld_H-1.5354_p1_branch2.res" u ($1):(-$3) w l ls 2 not, \
"stmfld_H-1.5354_p1_branch1.res" u ($1):(-$3) w l ls 1 not, \
"stmfld_H-1.5354_p1_branch2.res" u ($1):(-$3) w l ls 1 not, \
f(x) w l lt 1 lc 2 not

#"unstmfld_neg.res" w l lt 1 lc 1 not, \
#"stmfld.res" w l lt 1 lc 3 not, \
#"stmfld_neg.res" w l lt 1 lc 3 not, \
#"po.res" w p pt 3 pointsize 1 lc 2 not

set term pop
set out
