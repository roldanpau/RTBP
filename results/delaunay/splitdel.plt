set xlabel "e"
set ylabel "w_1"
#set out "splitdel.eps"
#set term post eps
plot [:0.89] "splittingdel.res" using 1:4 w l not
#set term pop
#set out
