#set title "Periodic orbit in Delaunay coordinates"
set xlabel "l"
set ylabel "L"
#set size square

set style arrow 1 nohead linetype 0 linewidth 1

# Plot horizontal line corresponding to trajectory of L1.
# This trajectory has L coordinate constant and equal to a^{1/2}, where 
# a is semi-major axis = radius of trajectory = x(L1) = 0.9324577513.
set arrow from 0,0.96563851999596619240 to 20,0.96563851999596619240 arrowstyle 1

#set out "trtbpdel.eps"
#set term post eps


#plot \

plot [0:20] [0.65:1] \
"trtbpdel_ell09.res" using 2:3 w p t "e=0.9", \
"trtbpdel_ell091.res" using 2:3 w p t "e=0.91", \
"trtbpdel_ell092.res" using 2:3 w p t "e=0.92"
#"trtbpdel_ell093.res" using 2:3 w l t "e=0.93", \
#"trtbpdel_ell094.res" using 2:3 w l t "e=0.94"

#"trtbpdel_ell01.res" using 2:3 w l t "e=0.1", \
#"trtbpdel_ell02.res" using 2:3 w l t "e=0.2", \
#"trtbpdel_ell03.res" using 2:3 w l t "e=0.3", \
#"trtbpdel_ell07.res" using 2:3 w l t "e=0.7", \
#"trtbpdel_ell08.res" using 2:3 w l t "e=0.8", \
#set term pop
#set out
