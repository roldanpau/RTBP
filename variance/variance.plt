#set out "variance.eps"
#set term postscript enhanced color eps

set hidden3d
set view 60, 350, 1, 1

set xlabel "I"
set ylabel "theta"
set zlabel offset -5,0 "s_0^2"

splot "variance.res" with lines palette not

set term pop
set out 

