set term tikz
set out "frtbpdel_dbg_zu_H-1.4854.tex"

set xlabel "$t$"
set ylabel "$dg/dt$"
#set title "Monitoring dg/dt for homoclinic orbit with H = -1.4854"

# Draw horizontal axis
#set arrow 1 from graph 0,first 0 to graph 1,first 0 nohead

set grid

plot "frtbpdel_dbg_zu_H-1.4854.err" u 1:3 w l not

set term pop
set out
