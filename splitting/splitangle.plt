unset label 

#set out "splitangle.eps"
#set term post eps
#set title "Splitting angle"
set xlabel "x"
set ylabel "px"
set label "z" at first -0.0870328, first 0

# Plot also homoclinic points?

plot [-0.08708:-0.08702] [-0.1:0.1] \
"unstmfld_ell085.res" w l lt 1 lc 1 not, \
"unstmfld_ell085.res" u 1:(-$2) w l lt 1 lc 2 not, \
"tanvec_u.res" u 1:2:($3)/25:($4)/25 w vectors lt 1 lc 1 t "wu", \
"tanvec_u.res" u 1:2:(-($3)/25):($4)/25 w vectors lt 1 lc 2 t "ws", \
"-" w p pt 3 pointsize 1 lc 3 not
-0.0870338 0
E

#set term pop
#set out
