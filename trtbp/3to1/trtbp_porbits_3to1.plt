set xlabel "x"
set ylabel "y"
#set title "Extremal periodic orbits in sideral coordinates"
set size square

#set out "trtbp_porbits_3to1.eps"
#set term post eps color

#plot \
#"trtbp_po_min.res" using ($2*cos($1)-$3*sin($1)):($2*sin($1)+$3*cos($1)) \
#w l t "H=-2.04", \
#"trtbp_po_max.res" using ($2*cos($1)-$3*sin($1)):($2*sin($1)+$3*cos($1)) \
#w l t "H=-1.57"

plot \
"trtbp_3to1.res" using 2:3 w l t "e=0", \
"trtbp_3to1_ell07.res" using 2:3 w l t "e=0.7", \
"L1.dat" not


#set term pop
#set out
