set xlabel "H"
set ylabel "T"
set y2label "L"
set ytics nomirror
set y2tics
set grid x y

set key t r spacing 3.0

#set out "porbits.eps" 
#set term post eps color
plot "porbits.res" using 1:(($2)-2.0*pi) w l lt 1 lc 1 t "T2XXXXXXX", \
"../Ldeviation/Ldeviation.res" w l lt 1 lc 2 t "L2XXXXXXX" axes x1y2

#set term pop
#set out 
