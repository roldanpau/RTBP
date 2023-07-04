set xlabel "H" 
unset ylabel

set key t r spacing 3.0

#set out "B_23.eps"
#set term postscript enhanced color eps

# plot Re(B^1), Im(B^1), Re(B^2), Im(B^2)

plot  \
"B_23.res" u 1:($5-$11) w l t "reB2-reB3", \
"B_23.res" u 1:($6-$12) w l t "imB2-imB3"


set term pop
set out
