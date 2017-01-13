set title "Upper branch of the manifolds"
set xlabel "H"
set ylabel "s"

set style line 1 lt 10 lc rgb "red" 
set style line 2 lt 0 lc rgb "blue" 
set style line 3 lt 3 lc rgb "green" 

f(x)=0.0

#set out "split_br1.eps"
#set term postscript eps color 
plot \
"splittings_SEC1_br1.res" w l ls 1 t "Channel 1", \
"splittings_SEC2_br1.res" w l ls 2 t "Channel 2", \
f(x) w l ls 3 not
#set term pop
#set out
