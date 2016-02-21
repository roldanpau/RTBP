unset logscale
set logscale y

#set title "Convergence between homoclinic and periodic trajectory"
set xlabel "N"	# multiples of the period
set ylabel "d"	# function dist(N)

set out "outer_circ_test.eps"
set term post eps colour

plot "outer_circ_test.res" w l lc palette z not

#plot "omega_test.res" every :::0::0 w l t "H=-1.6", \
#"omega_test.res" every :::1::1 w l t "H=-1.599", \
#"omega_test.res" every :::10::10 w l t "H=-1.586"

#"omega_test.res" every :::2::2 w l t "H=-1.598", \
#"omega_test.res" every :::3::3 w l t "H=-1.596", \
#"omega_test.res" every :::4::4 w l t "H=-1.594", \
#"omega_test.res" every :::5::5 w l t "H=-1.593"

set term pop
set out
