#set out "intersec.eps"
#set term post eps color
set xlabel "H"
set ylabel "x"

plot "intersec_p3.res" u 1:5 w l t "outer", \
"inner_sep/intersec_p3.res" u 1:5 w l t "inner"

#set term pop
#set out
