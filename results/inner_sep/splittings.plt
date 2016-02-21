#set out "splittings.eps"
#set term post eps color

set xlabel "H"
set ylabel "s"

set key l t

plot \
"../splitting.res" w l t "outer", \
"splitting.res" w l t "inner"
#set term pop
#set out
