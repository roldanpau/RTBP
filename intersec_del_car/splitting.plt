#set title "Upper branch of the manifolds"
set xlabel "H"
set ylabel "s"

f(x)=0.0

set out "splitting.eps"
set term postscript enhanced monochrome eps

set key left top horizontal

plot \
"splitting_SEC1_br1.res" w l t "s1p" smooth csplines, \
"splitting_SEC1_br2.res" w l t "s1m" smooth csplines, \
"splitting_SEC2_br1.res" w l t "s2p" smooth csplines, \
"splitting_SEC2_br2.res" w l t "s2m" smooth csplines, \
f(x) w l not

set term pop     
set out 
