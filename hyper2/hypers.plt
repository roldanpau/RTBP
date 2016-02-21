#set out "hypers.tex"
#set term epslatex

set xlabel "H"
#set logscale y
plot "hypers.res" u ($1):(log($2)) w l t "$\\lambda_u$", \
"hypers.res" u ($1):(log($3)) w l t "$\\lambda_s$"

#set term pop
#set out
