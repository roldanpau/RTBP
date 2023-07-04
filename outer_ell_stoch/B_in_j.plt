set xlabel "H" 
unset ylabel

set key t l spacing 3.0

#set out "B_in_j.eps"
#set term postscript enhanced monochrome eps

plot  \
"B_1.res" u 1:2 w l t "reBin", \
"B_1.res" u 1:3 w l t "imBin"


set term pop
set out
