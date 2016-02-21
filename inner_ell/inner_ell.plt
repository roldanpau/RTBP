#set out "inner_ell.eps"
#set term post enhanced color
set xlabel "H"
set ylabel "A^+(H)"
plot "inner_ell.res" u 1:5 w l t "|A^+|", \
 "inner_ell.res" u 1:4 w l t "Im(A^+)", \
 "inner_ell.res" u 1:3 w l t "Re(A^+)"
#set term pop
#set out
