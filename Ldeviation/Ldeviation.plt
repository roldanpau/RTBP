#set out "Ldeviation.eps"
#set term post eps

set xlabel "H"
set ylabel "Ldeviation"
plot [:-1.485] "Ldeviation.res" w l not

set term pop
set out
