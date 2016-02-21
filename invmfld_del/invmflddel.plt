unset label 

# fixed points
set label "p0" at first -1.042055571153243e+00, first 3.741457640450072e-01
set label "p1" at first 1.042055571153243e+00, first 3.741457640450072e-01

#set out "invmflddel_it_1.png"
#set term png
#set title "Invariant manifolds of fixed point"
set xlabel "g"
set ylabel "G"

set style arrow 1 nohead linetype 0 linewidth 1

# Plot vertical line corresponding to symmetry axis.
set arrow from 0,0 to 0,0.5 arrowstyle 1

plot [:] [:] \
"unstmfld.res" w l lt 1 lc 1 not, \
"unstmfld_neg.res" w l lt 1 lc 1 not, \
"stmfld.res" w l lt 1 lc 3 not, \
"stmfld_neg.res" w l lt 1 lc 3 not

pause -1

plot [-0.06:0.06] [0.34525:0.34535] \
"unstmfld.res" w l lt 1 lc 1 not, \
"unstmfld_neg.res" w l lt 1 lc 1 not, \
"stmfld.res" w l lt 1 lc 3 not, \
"stmfld_neg.res" w l lt 1 lc 3 not


#"stmfld_neg.res" w l lt 1 lc 3 not, \
#"po.res" w p pt 3 pointsize 1 lc 2 not
#set term pop
#set out
