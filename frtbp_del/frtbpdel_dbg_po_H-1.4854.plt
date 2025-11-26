#set out "frtbpdel_dbg_po_H-1.4854.eps"
#set term post eps color

set xlabel "t"
set ylabel "dg/dt"
#set title "Monitoring dg/dt for periodic orbit with H = -1.4854"

plot "frtbpdel_dbg_po_H-1.4854.err" u 1:3 w l not

set term pop
set out
