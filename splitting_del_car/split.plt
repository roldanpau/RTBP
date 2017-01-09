set xlabel "H"
set ylabel "s"

set style line 1 lt 10 lc rgb "red" 
set style line 2 lt 0 lc rgb "blue" 
set style line 3 lt 3 lc rgb "green" 

f(x)=0.0

#set out "split.eps"
#set term postscript eps color 
plot \
"splittings_unst_br1.res" w l ls 1 t "branch 1", \
"splittings_unst_br2.res" w l ls 2 t "branch 2", \
f(x) w l ls 3 not
#set term pop
#set out
