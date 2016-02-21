set xlabel "H"
set ylabel "s"

f(x)=0.0

#set out "split.eps"
#set term postscript eps color 
plot \
"splittings.res" w l t "s1", \
"splittings_branch2.res" w l t "s2", \
f(x) w l lt 3 lc 3 not
#set term pop
#set out
