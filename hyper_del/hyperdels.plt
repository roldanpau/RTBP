#set out "hyperdels.eps"
#set term post eps

set xlabel "e"
set ylabel "ln(lambda)"

plot "hyperdels.res" u ($1):(log($2)) w l not

#set term pop
#set out
