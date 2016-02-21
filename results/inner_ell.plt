#set out "inner_ell.eps"
#set term post eps 
set xlabel "H"
set ylabel "A"	# A_1^+(H)

set key b c spacing 3

# plot Re(A_1^+), Im(A_1^+)
plot \
 "inner_ell.res" u 1:2 w l t "reXXXXXX", \
 "inner_ell.res" u 1:3 w l t "imXXXXXX"
#set term pop
#set out
