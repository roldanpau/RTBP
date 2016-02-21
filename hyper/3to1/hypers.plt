#set out "hypers.eps"
#set term post eps

set xlabel "H"
set ylabel "lu"

#set logscale y
plot "hypers_3to1.res" u ($1):(log($2)) w l not

#set term pop
#set out
