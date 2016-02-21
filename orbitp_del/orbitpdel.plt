set xlabel "l"
set ylabel "L"

#set out "orbitpdel.eps"
#set term post eps color

# plot also $p_0$ (last), because orbitpdel does not print 0th iterate

plot [0:2*pi] \
"./unstmfld_it1.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_it1_neg.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_it1.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_it1_neg.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_it2.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_it2_neg.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_it2.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_it2_neg.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_it3.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_it3_neg.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_it3.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_it3_neg.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_it4.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_it4_neg.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_it4.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_it4_neg.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_it5.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_it5_neg.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_it5.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_it5_neg.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_it6.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_it6_neg.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_it6.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_it6_neg.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_it7.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_it7_neg.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_it7.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_it7_neg.res" w l lt 1 lw 3 lc 3 not, \
"./orbitpdel.res" w p pt 3 lc 2 not, \
"-" w p pt 3 lc 2 not 
0.000000000000000e+00 1.906394680726033e+00
E

#set term pop
#set out
