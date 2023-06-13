set xlabel "H" 
unset ylabel

set key t l spacing 3.0

#set out "B_12_difference.eps"
#set term postscript enhanced monochrome eps

# plot the difference Re(B^1) - Re(B^2), Im(B^1) - Im(B^2)

plot  \
"B_12.res" u 1:($2-$5) w l t "rediffXXXXXXXXXX", \
"B_12.res" u 1:($3-$6) w l t "imdiffXXXXXXXXXX"


set term pop
set out
