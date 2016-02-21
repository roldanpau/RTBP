#set out "Ldeviation.eps"
#set term post eps

set xlabel "H"
set ylabel "L"
plot "Ldeviation.res" w l not

#set term pop
#set out
