set term tikz
set out "B_out_j.tex"

set xlabel "$\\tJ$" 
unset ylabel

set key c r spacing 2.0

plot  \
"B_1.res" u 1:4 w l t "$\\Im{B_1^{\\out}}$", \
"B_2.res" u 1:4 w l t "$\\Im{B_2^{\\out}}$", \
"B_3.res" u 1:4 w l t "$\\Im{B_3^{\\out}}$", \
"B_4.res" u 1:4 w l t "$\\Im{B_4^{\\out}}$"

set term pop
set out
