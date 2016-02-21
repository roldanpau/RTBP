#set title "Periodic orbit in synodical coordinates"
set xlabel "x"
set ylabel "y"
#set size square
set label "z" at 1.399, 0.0
set label "P1" at 4.442, 0.0
set label "P2" at 6.113, -2
set label "P3" at 6.113, 0.0

set style arrow 1 nohead linetype 0 linewidth 1

#set out "trtbp.eps"
#set term post eps

# Plot line for label "P2"
set arrow from 5.620,0.0 to 6.19,-2 arrowstyle 1

plot "trtbp.res" using 2:3 w l not, \
"Piterates.res" w p pt 3 pointsize 1 not

#set term pop
#set out

#pause -1
#set out "trtbp_sideral.tex"
#set term epslatex color
#set title "Periodic orbit in sideral coordinates"
#set xlabel "X"
#set ylabel "Y"

# plot also trajectory of Kepler's problem (for comparison)
#plot \
#"trtbp.res" using ($2*cos($1)-$3*sin($1)):($2*sin($1)+$3*cos($1)) w l \
#   t "rtbp",\
#"ikepler.res" using ($2*cos($4)):($2*sin($4)) w l t "kepler"

#set term pop
#set out
