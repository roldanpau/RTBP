set term tikz
set out "omega_j.tex"

set xlabel "$\\tJ$"
set ylabel "$\\alpha_i$"

#set key left top
#set key spacing 3

set datafile missing "?"

plot \
"omega_neg.res" u 1:(-2.0*$2) t "$\\alpha_1$" w l, \
"omega_neg.res" u 1:(2.0*$4) t "$\\alpha_2$" w l, \
"omega_neg.res" u 1:(-2.0*$3) t "$\\alpha_3$" w l, \
"omega_neg.res" u 1:(2.0*$5) t "$\\alpha_4$" w l

set term pop
set out 
