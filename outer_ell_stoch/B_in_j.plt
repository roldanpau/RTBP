set term tikz
set out "B_in_j.tex"

set xlabel "$\\tJ$" 
unset ylabel

set key t l spacing 2.0 offset 5,0

plot  \
"B_1.res" u 1:2 w l t "$\\Re{B_1^\\inn}$", \
"B_1.res" u 1:3 w l t "$\\Im{B_1^\\inn}$"

#"B_2.res" u 1:2 w l t "$\\Re{B_2^\\inn}$", \
#"B_2.res" u 1:3 w l t "$\\Im{B_2^\\inn}$", \
#"B_3.res" u 1:2 w l t "$\\Re{B_3^\\inn}$", \
#"B_3.res" u 1:3 w l t "$\\Im{B_3^\\inn}$", \
#"B_4.res" u 1:2 w l t "$\\Re{B_4^\\inn}$", \
#"B_4.res" u 1:3 w l t "$\\Im{B_4^\\inn}$"


set term pop
set out
