set xlabel "H" 
unset ylabel

set key t r spacing 3.0

set out "B_12.eps"
set term postscript enhanced color eps

# plot Re(B^1), Im(B^1), Re(B^2), Im(B^2)

#plot [-1.81:-1.56] \
#"./inner_sep/B_f.res" u 1:2 w l t "reBfXXXXXX", \
#"./inner_sep/B_f.res" u 1:3 w l t "imBfXXXXXX", \
#"./B_b.res" u 1:2 w l t "reBbXXXXXX", \
#"./B_b.res" u 1:3 w l t "imBbXXXXXX"

plot  \
"B_12.res" u 1:2 w l t "reB1", \
"B_12.res" u 1:3 w l t "imB1", \
"B_12.res" u 1:5 w l t "reB2", \
"B_12.res" u 1:6 w l t "imB2", \
"B_12.res" u 1:($2-$5) w l t "reB1-reB2", \
"B_12.res" u 1:($3-$6) w l t "imB1-imB2"


set term pop
set out
