# Check how far is period of periodic orbits from 2\pi

set xlabel "H"
set ylabel "T"
set y2label "L"
set ytics nomirror
set y2tics
set grid x y

set key t l spacing 3.0

#set out "porbits_3to1.eps" 
#set term post eps color
plot "porbits_3to1.res" using 1:(($2)-2.0*pi) w l lt 1 lc 1 t "T2XXXXXXX"

#"../Ldeviation/Ldeviation.res" w l lt 1 lc 2 t "L2XXXXXXX" axes x1y2

#set term pop
#set out 
