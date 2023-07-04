set xlabel "H" 
unset ylabel

set key t l spacing 3.0

#set out "B_out_j.eps"
#set term postscript enhanced monochrome eps

plot  \
"B_1.res" u 1:4 w l t "imBout1", \
"B_2.res" u 1:4 w l t "imBout2", \
"B_3.res" u 1:4 w l t "imBout3", \
"B_4.res" u 1:4 w l t "imBout4"

set term pop
set out
