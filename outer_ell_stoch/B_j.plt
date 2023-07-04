set xlabel "H" 
unset ylabel

set key t r spacing 3.0

#set out "B_j.eps"
#set term postscript enhanced color eps

# plot Re(B^1), Im(B^1), Re(B^2), Im(B^2)

plot  \
"B_1.res" u 1:6 w l t "imB1", \
"B_2.res" u 1:6 w l t "imB2", \
"B_3.res" u 1:6 w l t "imB3", \
"B_4.res" u 1:6 w l t "imB4"

#"B_12.res" u 1:($5-$11) w l t "reB1-reB2", \
#"B_12.res" u 1:($6-$12) w l t "imB1-imB2"


set term pop
set out
