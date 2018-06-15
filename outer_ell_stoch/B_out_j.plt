set xlabel "H" 
unset ylabel

set key t l spacing 3.0

#set out "B_out_j.eps"
#set term postscript enhanced monochrome eps

plot  \
"B_out_1.res" u 1:2 w l t "reB1p", \
"B_out_1.res" u 1:3 w l t "imB1p", \
"B_out_2.res" u 1:2 w l t "reB2p", \
"B_out_2.res" u 1:3 w l t "imB2p"


set term pop
set out
