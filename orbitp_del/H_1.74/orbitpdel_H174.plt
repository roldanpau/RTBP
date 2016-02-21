set xlabel "l"
set ylabel "L"

set label "p0" at first 3.151257033096044e-12, first 1.909711119665817e+00
set label "p1" at first 8.972360712083836e-01, first 1.912886782966373e+00
set label "p2" at first 1.794995867222006e+00, first 1.912862153280634e+00
set label "p3" at first 2.692727947084270e+00, first 1.912858911679164e+00
set label "p4" at first 3.590457360108106e+00, first 1.912858911679167e+00
set label "p5" at first 4.488189439964807e+00, first 1.912862153280978e+00
set label "p6" at first 5.385949235977967e+00, first 1.912886782966720e+00

# plot vertical line corresponding to g=pi axis of symmetry

set arrow from graph 0.5,0 to graph 0.5,1 nohead 

set label "z1" at first pi, first 1.922
set label "z2" at first pi, first 1.904

#set out "orbitpdel_H174.eps"
#set term post eps color

plot [0:2*pi] \
"./unstmfld_p0.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_neg_p0.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_p0.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_neg_p0.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_p1.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_neg_p1.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_p1.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_neg_p1.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_p2.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_neg_p2.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_p2.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_neg_p2.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_p3.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_neg_p3.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_p3.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_neg_p3.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_p4.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_neg_p4.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_p4.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_neg_p4.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_p5.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_neg_p5.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_p5.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_neg_p5.res" w l lt 1 lw 3 lc 3 not, \
"./unstmfld_p6.res" w l lt 1 lw 3 lc 1 not, \
"./unstmfld_neg_p6.res" w l lt 1 lw 3 lc 1 not, \
"./stmfld_p6.res" w l lt 1 lw 3 lc 3 not, \
"./stmfld_neg_p6.res" w l lt 1 lw 3 lc 3 not, \
"./po.res" w p pt 3 lc 2 not, \
"-" w p pt 3 lc 2 not 
3.1416 1.92251
3.1416 1.90343
E

#set term pop
#set out
