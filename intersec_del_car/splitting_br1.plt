#set title "Upper branch of the manifolds"
set xlabel "H"
set ylabel "s"

#set style line 1 lt 10 lc rgb "red" 
#set style line 2 lt 0 lc rgb "blue" 
#set style line 3 lt 3 lc rgb "green" 

set style line 1 lt 10
set style line 2 lt 0
set style line 3 lt 10

f(x)=0.0

#set out "split_br1.eps"
#set term postscript eps 

plot \
"splitting_SEC1_br1.res" w l ls 1 not smooth csplines, \
"splitting_SEC2_br1.res" w l ls 2 not smooth csplines, \
f(x) w l ls 3 not

#set term pop     
#set out 
