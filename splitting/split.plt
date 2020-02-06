set xlabel "H"
set ylabel "s"

f(x)=0.0

#set out "split.eps"
#set term postscript eps color 
plot \
"splittings_unst_br1.res" w l t "su1", \
"splittings_unst_br2.res" w l t "su2", \
"splittings_branch2.res" w l t "s2", \
f(x) w l lt 3 lc 3 not
set term pop
set out
