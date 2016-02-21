set out "diffs.eps"
set term post eps

set xlabel "H"
set ylabel "srad"

plot [-1.81:-1.8] "diffs.res" u 1:2 w l t "s", \
"diffs.res" u 1:($2-$3) w l t "err"
set term pop
set out
