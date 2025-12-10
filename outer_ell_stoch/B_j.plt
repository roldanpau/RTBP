set term tikz
set out "B_j.tex"

set xlabel "$\\tJ$" 
unset ylabel

set key c r spacing 2.0

plot  \
"B_1.res" u 1:6 w l t "$\\Im{B_1}$", \
"B_2.res" u 1:6 w l t "$\\Im{B_2}$", \
"B_3.res" u 1:6 w l t "$\\Im{B_3}$", \
"B_4.res" u 1:6 w l t "$\\Im{B_4}$"

#"B_12.res" u 1:($5-$11) w l t "reB1-reB2", \
#"B_12.res" u 1:($6-$12) w l t "imB1-imB2"


set term pop
set out
