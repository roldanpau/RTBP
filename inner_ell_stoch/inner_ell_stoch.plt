#set term postscript enhanced monochrome eps
#set out "inner_ell_stoch.eps"

set linestyle 1 linetype 1
set linestyle 2 linetype 2

set xlabel "H"
set ylabel "A(H)"
plot \
 "inner_ell_stoch.res" u 1:2 w l ls 1 t "Re(A)", \
 "inner_ell_stoch.res" u 1:3 w l ls 2 t "Im(A)"
set term pop
set out
