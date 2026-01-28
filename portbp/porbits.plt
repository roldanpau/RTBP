#set term tikz
#set out "porbits.tex"

set xlabel "$\\tJ$"
set ylabel "$T$"
set y2label "$L$"
set ytics nomirror
set y2tics
set grid x y

set key t r spacing 3.0

plot [:-1.485] "porbits.res" using 1:(($2)-2.0*pi) w l \
t "$|T_\\tJ - 2\\pi|$", \
"../Ldeviation/Ldeviation.res" w l \
t "$|\\overline{L_\\tJ} - 3^{-1/3}|$" axes x1y2

#plot [:-1.485] "porbits.res" using 1:(($2)-2.0*pi) w l \
#t "$|T_\\tJ - 2\\pi|$", \
#"../Ldeviation/Ldeviation.res" using 1:3 w l \
#t "$|\\overline{a_\\tJ} - 3^{-2/3}|$" axes x1y2

#set term pop
#set out 
