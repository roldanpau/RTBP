set xlabel "$H$"
set ylabel "$\\mu T_0$"

#set out "inner_circ.tex" 
#set term epslatex 
plot "inner_circ.res" using 1:(($2)-14.0*pi) w l not
#set term pop
#set out 
