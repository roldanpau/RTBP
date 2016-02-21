# Symmetric invariant manifolds

#set out "invmfld_sym.tex"
#set term epslatex color
set title "Symmetric invariant manifolds of fixed point"
set xlabel "$x$"
set ylabel "$p_x$"

plot "unstmfld.res" w l, "stmfldup.res" w l
#set term pop
#set out
pause -1

#set out "invmfld_zoom.tex"
#set term epslatex color
#set title "Invariant manifolds of fixed point (zoom)"
#set xlabel "$x$"
#set ylabel "$p_x$"
#
#plot [2.19:3.10] [:-0.38] "unstmfld.res" w l, "stmfld.res" w l
#set term pop
#set out
#pause -1
