set xlabel "H"
set ylabel "s"

f(x)=0.0

#set out "split.eps"
#set term postscript eps color 
plot \
"splittings_unst_br1.res" w l t "u1", \
"splittings_unst_br2.res" w l t "u2", \
"splittings_st_br1.res" w l t "s1", \
"splittings_st_br2.res" w l t "s2", \
f(x) w l lt 3 lc 3 not
set term pop
set out
